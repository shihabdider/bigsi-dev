/* Queries a bigsi matrix of buckets from a reference sequence
 *
 * Input: A query sequence as a string
 * Output: An object containing hits per bucket where hits are either
 *  a) number of query fragments found in bucket
 *  b) number of minimizers found in bucket
 *
 *  A query fragment represents an exact match of length equal to the fragment 
 *  (e.g 2.5Kbp). 
 *
 *  A minimizer represents an exact match of (w + k - 1) bp, i.e 
 *  an exact match of (w + k - 1) bps guarantees a shared minimizer with a 
 *  false positive rate equal to that of the Bloom Filter (generally too small 
 *  to be of concern). w, indicates the window size used to construct the 
 *  minimizer and is set 1000, meaning each minimizer "represents" 500 k-mers.
 *  
 */

const cdf = require('binomial-cdf')
const matrix = require('matrix-js')
const BitSet = require('bitset')
const fs = require('fs')
const utils = require('./utils.js')
const quantile = require( '@stdlib/stats-base-dists-binomial-quantile' );

// One function is used for fragmenting and winnowing to prevent double 
// iteration
async function winnowQueryFragments(querySeq, fragmentSize=5000){

    const querySize = querySeq.length
    const queryWindowSize = 100

    const queryFragmentsMinimizers = []

    if (fragmentSize == 0){ //for non-frag queries
        const queryMinimizers = utils.extractMinimizers(querySeq, queryWindowSize)
        queryFragmentsMinimizers.push(queryMinimizers)
    } else {

        if (querySize <= 300_000 && querySize >= fragmentSize){

            for (let n=0; n < Math.floor(querySize/fragmentSize); n++){
                const start = n*fragmentSize
                const end = Math.min(start + fragmentSize, querySize)
                const queryFragmentSeq = querySeq.slice(start, end)
                const queryFragmentMinimizers = utils.extractMinimizers(queryFragmentSeq, 1000)
                queryFragmentsMinimizers.push(queryFragmentMinimizers)
            }
        } else { console.error("Inappropriate query size") }
    }

    return queryFragmentsMinimizers
}

function makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize){

    const fragmentsBloomFilters = []
    for (let fragmentMinimizerSet of queryFragmentsMinimizers){
        const bf = utils.makeMinimizersBloomFilter(fragmentMinimizerSet, bloomFilterSize)
        fragmentsBloomFilters.push(bf)
    }

    return fragmentsBloomFilters
}

function getBloomFilterSetBitsIndices(queryBF){
    return queryBF.reduce((indices, number, index) => {
        if (number != 0) indices.push(index)
        return indices
    }, [])
}

/**
 * @param { Uint16Array } bigsi - flattened bigsi 
 * @param { array of numbers } rowFilter - array of row numbers for query rows
 * @param { number } numCols - number of columns in bigsi
 *
 * @returns { matrix } submatrix - query rows submatrix
 */
function getBinaryBigsiSubmatrix(bigsi, rowFilter, numCols){
    const submatrixRows = []

    const intSize = 16
    const numInts = numCols/intSize

    for (const rowNum of rowFilter){

        const offsetStart = rowNum*numInts
        const offsetEnd = offsetStart + numInts

        const rowInts = Array.from(bigsi.subarray(offsetStart, offsetEnd))
        const rowBitStrings = rowInts.map((num) => utils.zeroPadBitstring(num.toString(2), 16))
        const rowBitString = rowBitStrings.join('')

        const row = rowBitString.split('').map(Number)
        submatrixRows.push(row)
    }

    const submatrix = matrix(submatrixRows)

    return submatrix
}

function computeSubmatrixHits(submatrix, bigsiHits, numBuckets) {
    const numRows = submatrix.size()[0]
    const submatrix_T = submatrix.trans()
    const hammingWeights = []
    for (rowNum=0; rowNum < numBuckets; rowNum++){
        const row = submatrix_T[rowNum]
        const bs = new BitSet(row.join(''))
        const weight = bs.cardinality()
        if (weight == numRows) {
            //console.log('row num', weight)
            if (rowNum in bigsiHits){
                bigsiHits[rowNum]['hits'] += 1
            } else {
                bigsiHits[rowNum] = {'hits': 1}
            }
        }
    }
}

function computeLowerBoundContainmentScore(containmentScore, 
    numMinimizersInQuery, confidenceInterval) {
    // begin search from x = s * containment score
    let x = quantile(confidenceInterval, numMinimizersInQuery, containmentScore)

    const lowerBoundContainmentScore = Math.min(x / numMinimizersInQuery, 1);
    return lowerBoundContainmentScore; 
}

// Containment score is the Jaccard containment identity:
// Hamming weight of submatrix columns divided by
// number of minimizers inserted into query Bloom Filter
function computeQueryContainmentScores(submatrix, bigsiHits, bloomFilterSize, subrate) {
    const kmerLength = 16
    const queryNumBitsSet = submatrix.size()[0]
    const submatrix_T = submatrix.trans()
    const hammingWeights = []
    for (const row of submatrix_T){
        const bs = new BitSet(row.join(''))
        const weight = bs.cardinality()
        hammingWeights.push(weight)
    }

    for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
        let numIntersections = hammingWeights[bucketNum]
        if (numIntersections > 0) {
            const containmentScore = numIntersections/queryNumBitsSet
            const lowerContainment = computeLowerBoundContainmentScore(containmentScore, queryNumBitsSet, 0.9999999)
            const errorRate = Math.max(-1/kmerLength * Math.log(lowerContainment), 0)
            if (errorRate <= subrate) {
                //console.log(hammingWeights[bucketNum])
                const percentMatch = 100*(1 - errorRate)
                bigsiHits[bucketNum] = {'percent match': percentMatch}
            }
        }
    }
}


/** 
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { string } numCols - number of columns in bigsi
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryBinaryBigsi(bigsiArray, queryFragmentsBloomFilters, numCols, bloomFilterSize, subrate){

    const bigsiHits = {}

    const numFragments = queryFragmentsBloomFilters.length
    //console.log('number of query fragments: ', numFragments)

    for (const bloomFilter of queryFragmentsBloomFilters){
        const queryBFSetBitsIndices = getBloomFilterSetBitsIndices(bloomFilter)
        
        const querySubmatrix = await getBinaryBigsiSubmatrix(bigsiArray, queryBFSetBitsIndices, numCols)

        if (numFragments == 1){
            computeQueryContainmentScores(querySubmatrix, bigsiHits, bloomFilterSize, subrate)
        } else {
            computeSubmatrixHits(querySubmatrix, bigsiHits, numCols)
        }
    }

    for (const bucketId in bigsiHits) {
        bigsiHits[bucketId]['score'] = `${bigsiHits[bucketId]['hits']}/${numFragments}`;
    }

    return bigsiHits
}

async function main(querySeq, bigsiPath, bigsiConfigPath, subrate) {
    const bigsiBuffer = fs.readFileSync(bigsiPath)
    let bigsiArray = new Uint16Array(bigsiBuffer.buffer, bigsiBuffer.byteOffset, bigsiBuffer.length / 2);
    //console.log('bigsiArray size: ', bigsiArray.length)

    const bigsiDims = require(bigsiConfigPath)
    const numCols = bigsiDims.cols
    const bloomFilterSize = bigsiDims.rows
    //const queryFragmentsMinimizers = await winnowQueryFragments(querySeq)
    // Test: non-frag query
    const isQuerySeqRightSize = querySeq.length >= 1000 && querySeq.length <= 300000
    if (true) {
        const fragmentSizeZero = 0
        const queryFragmentsMinimizers = await winnowQueryFragments(querySeq, fragmentSizeZero)
        const queryFragmentsBloomFilters = await makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)
        //const queryMinimizers = queryFragmentsMinimizers[0]
        //const numHashes = 5
        //const queryRowFilters = await utils.makeQueryRowFilters(queryMinimizers, bloomFilterSize, numHashes)

        const filteredBigsiHits = await queryBinaryBigsi(bigsiArray, queryFragmentsBloomFilters, numCols, bloomFilterSize, subrate)

        return filteredBigsiHits
    } else {
        console.log('Query must be of length between 5Kb and 300Kb')
    }
}


module.exports = {
    winnowQueryFragments: winnowQueryFragments,
    makeFragmentsBloomFilters: makeFragmentsBloomFilters,
    getBloomFilterSetBitsIndices: getBloomFilterSetBitsIndices,
    getBinaryBigsiSubmatrix: getBinaryBigsiSubmatrix,
    computeSubmatrixHits: computeSubmatrixHits,
    computeQueryContainmentScores: computeQueryContainmentScores, 
    queryBinaryBigsi: queryBinaryBigsi,
    main: main
}
