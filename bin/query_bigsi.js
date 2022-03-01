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
const helper = require('./helper.js')

// One function is used for fragmenting and winnowing to prevent double 
// iteration
async function winnowQueryFragments(querySeq, fragmentSize=5000){

    const querySize = querySeq.length
    const queryWindowSize = 100

    const queryFragmentsMinimizers = []

    if (fragmentSize == 0){ //for non-frag queries
        const queryMinimizers = helper.extractMinimizers(querySeq, queryWindowSize)
        queryFragmentsMinimizers.push(queryMinimizers)
    } else {

        if (querySize <= 300_000 && querySize >= fragmentSize){

            for (let n=0; n < Math.floor(querySize/fragmentSize); n++){
                const start = n*fragmentSize
                const end = Math.min(start + fragmentSize, querySize)
                const queryFragmentSeq = querySeq.slice(start, end)
                const queryFragmentMinimizers = helper.extractMinimizers(queryFragmentSeq, 1000)
                queryFragmentsMinimizers.push(queryFragmentMinimizers)
            }
        } else { console.error("Inappropriate query size") }
    }

    return queryFragmentsMinimizers
}

function makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize){

    const fragmentsBloomFilters = []
    for (let fragmentMinimizerSet of queryFragmentsMinimizers){
        const bf = helper.makeMinimizersBloomFilter(fragmentMinimizerSet, bloomFilterSize)
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

    const numSeqs = numCols/16

    for (const rowNum of rowFilter){

        const offsetStart = rowNum*numSeqs
        const offsetEnd = offsetStart + numSeqs

        const rowInts = Array.from(bigsi.subarray(offsetStart, offsetEnd))
        const rowBitStrings = rowInts.map((num) => helper.zeroPadBitstring(num.toString(2), 16))
        const rowBitString = rowBitStrings.join('')

        // Front padding ensures all columns are accounted for
        const row = rowBitString.split('').map(Number)
        submatrixRows.push(row)
    }

    const submatrix = matrix(submatrixRows)

    return submatrix
}

// hexBigsi is an array of hexstrings
function getHexBigsiSubmatrix(hexBigsi, rowFilter){
    const submatrixRows = []

    for (const rowNum of rowFilter){

        const rowHex = hexBigsi[rowNum]
        const rowBitString = parseInt(rowHex, 16).toString(2).padStart(16, '0')
        // Front padding ensures all columns are accounted for
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
            console.log('row num', weight)
            if (rowNum in bigsiHits){
                bigsiHits[rowNum]['hits'] += 1
            } else {
                bigsiHits[rowNum] = {'hits': 1}
            }
        }
    }
}

/**
 * @param { matrix } submatrix - array of submatrix for query rows 
 * of bigsi
 * @param { bigsiHits } bigsiHits - empty object for storing returned hits 
 *
 * @returns - bigsiHits with updated hits attributes 
 */
//function computeSubmatrixHits(submatrix, bigsiHits, numBuckets){
//    let bucketHits = (new BitSet).flip() // init bitset to all 1s for bitwise AND
//    const numRows = submatrix.size()[0]
//    for (let row = 0; row < numRows; row++) {
//        const rowBitString = submatrix(row).join('')
//        const rowBS = new BitSet(rowBitString)
//        bucketHits = bucketHits.and(rowBS)
//    }
//
//    // Use bigsiHits to compute the total number of buckets instead of bitset 
//    // length because the latter truncates leading zeros when converted to 
//    // a string
//    //const numBuckets = Object.keys(bigsiHits).length
//    
//    // bitset.toArray returns an array of indices (i) corresponding to the set bits
//    // to convert this to bucket numbers (B - 1) - i, where B is 
//    // the total number of buckets
//    const hitsBucketNums = bucketHits.toArray().map( value => (numBuckets - 1) - value )
//
//    for (let i = 0; i < hitsBucketNums.length; i++){
//        const bucketNum = hitsBucketNums[i]
//        if (bucketNum in bigsiHits){
//            bigsiHits[bucketNum]['hits'] += 1
//        } else {
//            bigsiHits[bucketNum] = {'hits': 1}
//        }
//
//    }
//}

// Containment score is the Jaccard containment identity:
// Hamming weight of submatrix columns divided by
// number of minimizers inserted into query Bloom Filter
// (using union approximation for multiple hashes)
//function computeQueryContainmentScores(queryBF, bigsi, bigsiHits, bloomFilterSize) {
//    const queryMinimizerCount = submatrix.size()[0]
//    const submatrix_T = submatrix.trans()
//    const hammingWeights = []
//    for (const row of submatrix_T){
//        const bs = new BitSet(row.join(''))
//        const weight = bs.cardinality()
//        hammingWeights.push(weight)
//    }
//
//    //const numHashes = 1
//    const numHashes = 5
//    for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
//        let numIntersections = hammingWeights[bucketNum]
//        const numQueryMinimizers = (-1*bloomFilterSize/numHashes)*Math.log(1 - queryMinimizerCount/bloomFilterSize)
//        const numBucketMinimizers = (-1*bloomFilterSize/numHashes)*Math.log(1 - bucketBFlength/bloomFilterSize)
//        const numUnionMinimizers = (-1*bloomFilterSize/numHashes)*Math.log(1 - numUnionBits/bloomFilterSize)
//        if (numHashes > 1) {
//            numIntersections = numQueryMinimizers + numBucketMinimizers - numUnionMinimizers
//        }
//        const containmentScore = numIntersections/numQueryMinimizers
//        const errorRate = -1/16 * Math.log(containmentScore)
//        if (errorRate <= 0.2) {
//            bigsiHits[bucketNum] = {'error rate': errorRate}
//        }
//    }
//    //let maxContainment = {'bucketNum': 0, 'containment': 0}
//    //for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
//    //    const containmentScore = hammingWeights[bucketNum]/queryMinimizerCount
//    //    if (containmentScore > maxContainment.containment) {
//    //        maxContainment.bucketNum = bucketNum
//    //        maxContainment.containment = containmentScore
//    //    }
//    //}
//    //bigsiHits[maxContainment.bucketNum] = {'containment': maxContainment.containment}
//
//}

// Containment score is the Jaccard containment identity:
// Hamming weight of submatrix columns divided by
// number of minimizers inserted into query Bloom Filter
// (using intersection approximation for multiple hashes)
function computeQueryContainmentScores(submatrix, bigsiHits, bloomFilterSize) {
    const queryNumBitsSet = submatrix.size()[0]
    const submatrix_T = submatrix.trans()
    const hammingWeights = []
    for (const row of submatrix_T){
        const bs = new BitSet(row.join(''))
        const weight = bs.cardinality()
        hammingWeights.push(weight)
    }

    const numHashes = 1
    //const numHashes = 5
    //const numQueryMinimizers = (-1*bloomFilterSize/numHashes)*Math.log(1 - queryNumBitsSet/bloomFilterSize)
    const numQueryMinimizers = queryNumBitsSet
    for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
        let numIntersections = hammingWeights[bucketNum]
        if (numHashes > 1) {
            numIntersections = (-1*queryNumBitsSet/numHashes)*Math.log(1 - hammingWeights[bucketNum]/queryNumBitsSet)
        }
        const containmentScore = numIntersections/numQueryMinimizers
        const errorRate = -1/16 * Math.log(containmentScore)
        if (errorRate <= 0.05) {
            console.log(hammingWeights[bucketNum])
            bigsiHits[bucketNum] = {'error rate': errorRate}
        }
    }
    //let maxContainment = {'bucketNum': 0, 'containment': 0}
    //for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
    //    const containmentScore = hammingWeights[bucketNum]/queryMinimizerCount
    //    if (containmentScore > maxContainment.containment) {
    //        maxContainment.bucketNum = bucketNum
    //        maxContainment.containment = containmentScore
    //    }
    //}
    //bigsiHits[maxContainment.bucketNum] = {'containment': maxContainment.containment}

}

function queryHexBigsi(hexBigsi, queryFragmentsBloomFilters, bloomFilterSize){
    const bigsiHits = {}

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query fragments: ', numFragments)

    for (const bloomFilter of queryFragmentsBloomFilters){
        const queryBFSetBitsIndices = getBloomFilterSetBitsIndices(bloomFilter._filter)
        
        const querySubmatrix = getHexBigsiSubmatrix(hexBigsi, queryBFSetBitsIndices)

        if (numFragments == 1){
            computeQueryContainmentScores(querySubmatrix, bigsiHits, bloomFilterSize)
        } else {
            computeSubmatrixHits(querySubmatrix, bigsiHits, numCols)
        }
    }

    for (const bucketId in bigsiHits) {
        bigsiHits[bucketId]['score'] = `${bigsiHits[bucketId]['hits']}/${numFragments}`;
    }

    return bigsiHits
}

/** 
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { string } numCols - number of columns in bigsi
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryBinaryBigsi(bigsiArray, queryFragmentsBloomFilters, numCols, bloomFilterSize){

    const bigsiHits = {}

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query fragments: ', numFragments)

    for (const bloomFilter of queryFragmentsBloomFilters){
        const queryBFSetBitsIndices = getBloomFilterSetBitsIndices(bloomFilter._filter)
        
        const querySubmatrix = await getBinaryBigsiSubmatrix(bigsiArray, queryBFSetBitsIndices, numCols)

        if (numFragments == 1){
            computeQueryContainmentScores(querySubmatrix, bigsiHits, bloomFilterSize)
        } else {
            computeSubmatrixHits(querySubmatrix, bigsiHits, numCols)
        }
    }

    for (const bucketId in bigsiHits) {
        bigsiHits[bucketId]['score'] = `${bigsiHits[bucketId]['hits']}/${numFragments}`;
    }

    return bigsiHits
}

/** 
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { string } numCols - number of columns in bigsi
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
//async function queryBinaryBigsi(bigsiArray, queryRowFilters, numCols, bloomFilterSize){
//
//    const bigsiHits = {}
//
//    for (const queryBFSetBitsIndices of queryRowFilters){
//        const querySubmatrix = await getBinaryBigsiSubmatrix(bigsiArray, queryBFSetBitsIndices, numCols)
//        computeSubmatrixHits(querySubmatrix, bigsiHits, numCols)
//    }
//
//    for (const bucketId in bigsiHits) {
//        bigsiHits[bucketId]['score'] = bigsiHits[bucketId]['hits']/queryRowFilters.length
//    }
//
//    return bigsiHits
//}

async function main(querySeq, bigsiPath, bigsiConfigPath) {
    const bigsiBuffer = fs.readFileSync(bigsiPath)
    let bigsiArray = new Uint16Array(bigsiBuffer.buffer, bigsiBuffer.byteOffset, bigsiBuffer.length / 2);
    console.log('bigsiArray size: ', bigsiArray.length)

    const bigsiDims = require(bigsiConfigPath)
    const numCols = bigsiDims.cols
    const bloomFilterSize = bigsiDims.rows
    //const queryFragmentsMinimizers = await winnowQueryFragments(querySeq)
    // Test: non-frag query
    const fragmentSizeZero = 0
    const queryFragmentsMinimizers = await winnowQueryFragments(querySeq, fragmentSizeZero)
    const queryFragmentsBloomFilters = await makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)
    //const queryMinimizers = queryFragmentsMinimizers[0]
    //const numHashes = 5
    //const queryRowFilters = await helper.makeQueryRowFilters(queryMinimizers, bloomFilterSize, numHashes)

    const filteredBigsiHits = await queryBinaryBigsi(bigsiArray, queryFragmentsBloomFilters, numCols, bloomFilterSize)

    return filteredBigsiHits
}

async function run() {
    if (require.main === module) {
        const { argv } = require('yargs')
            .scriptName('query_bigsi')
            .usage('Usage: $0 --query (path to query fasta) --bigsi (path to BIGSI file) --config (path to bigsi config)')
            .option('query', {
                alias: 'q',
                describe: 'Path to fasta file of query sequence',
                demandOption: 'Fasta file is required',
                type: 'string',
                nargs: 1,
            })
            .option('bigsi', {
                alias: 'b',
                describe: 'Path to BIGSI file',
                demandOption: 'BIGSI file is required',
                type: 'string',
                nargs: 1,
            })
            .option('config', {
                alias: 'c',
                describe: 'Path to BIGSI config file',
                demandOption: 'BIGSI config file is required',
                type: 'string',
                nargs: 1,
            })

        const fai = `${argv.query}.fai`
        const query = await helper.loadFasta(argv.query, fai)
        //console.log('Query sequence loaded...', query)
        const seqNames = await query.getSequenceList()
        console.log('Sequence name: ', seqNames)
        const querySeq = await query.getSequence(seqNames[0])
        const bigsiPath = argv.bigsi
        const bigsiConfigPath = argv.config
        const hits = await main(querySeq, bigsiPath, bigsiConfigPath);
        console.log(hits)
    }
}

run()

module.exports = {
    winnowQueryFragments: winnowQueryFragments,
    makeFragmentsBloomFilters: makeFragmentsBloomFilters,
    getBloomFilterSetBitsIndices: getBloomFilterSetBitsIndices,
    getBinaryBigsiSubmatrix: getBinaryBigsiSubmatrix,
    getHexBigsiSubmatrix: getHexBigsiSubmatrix,
    computeSubmatrixHits: computeSubmatrixHits,
    computeQueryContainmentScores: computeQueryContainmentScores, 
    queryBinaryBigsi: queryBinaryBigsi,
    queryHexBigsi: queryHexBigsi,
    main: main
}
