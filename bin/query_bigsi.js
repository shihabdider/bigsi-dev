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

const matrix = require('matrix-js')
const BitSet = require('bitset')
const fs = require('fs')
const utils = require('./utils.js')
const quantile = require( '@stdlib/stats-base-dists-binomial-quantile' );

function getBloomFilterSetBitsIndices(bloomFilter){
    return bloomFilter.reduce((indices, number, index) => {
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


function computeLowerBoundContainmentScore(containmentScore, 
    numMinimizersInQuery, confidenceInterval) {
    let x = quantile(confidenceInterval, numMinimizersInQuery, containmentScore)

    const lowerBoundContainmentScore = Math.min(x / numMinimizersInQuery, 1);
    return lowerBoundContainmentScore; 
}

// Containment score is the Jaccard containment identity:
// Hamming weight of submatrix columns divided by
// number of minimizers inserted into query Bloom Filter
function computeQueryContainmentScores(submatrix, bigsiHits, bloomFilterSize, subrate) {
    const kmerLength = 16
    const numMinimizersInQuery = submatrix.size()[0]
    const submatrix_T = submatrix.trans()
    const hammingWeights = []
    for (const row of submatrix_T){
        const bs = new BitSet(row.join(''))
        const weight = bs.cardinality()
        hammingWeights.push(weight)
    }

    for (let bucketNum = 0; bucketNum < hammingWeights.length; bucketNum++){
        let numMatchingMinimizers = hammingWeights[bucketNum]
        if (numMatchingMinimizers > 0) {
            let containmentScore = numMatchingMinimizers/numMinimizersInQuery
            const containmentBias = -0.1419
            containmentScore += containmentBias
            const lowerContainment = computeLowerBoundContainmentScore(containmentScore, numMinimizersInQuery, 0.99)
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
 * @param { Array } queryMask - bloom filter of the query sequence as an array
 * @param { string } numBins - number of bins in bigsi (corresponding to number 
 * of columns)
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryBinaryBigsi(bigsiArray, queryMask, numBins, subrate){

    const bigsiHits = {}

    const queryMaskSize = queryMask.length
    const queryBFSetBitsIndices = getBloomFilterSetBitsIndices(queryMask)
    const querySubmatrix = await getBinaryBigsiSubmatrix(bigsiArray, queryBFSetBitsIndices, numBins)

    computeQueryContainmentScores(querySubmatrix, bigsiHits, queryMaskSize, subrate)

    return bigsiHits
}

async function main(querySeq, bigsiPath, bigsiConfigPath, subrate) {
    const bigsiBuffer = fs.readFileSync(bigsiPath)
    let bigsiArray = new Uint16Array(bigsiBuffer.buffer, bigsiBuffer.byteOffset, bigsiBuffer.length / 2);
    //console.log('bigsiArray size: ', bigsiArray.length)

    const bigsiDims = require(bigsiConfigPath)
    const numBins = bigsiDims.cols
    const bloomFilterSize = bigsiDims.rows
    const queryWindowSize = 100

    const isQuerySeqRightSize = querySeq.length >= 1000 && querySeq.length <= 300000
    if (isQuerySeqRightSize) {
        const queryMinimizers = utils.extractMinimizers(querySeq, queryWindowSize)
        const queryMask = utils.makeMinimizersBloomFilter(queryMinimizers, bloomFilterSize)
        const filteredBigsiHits = await queryBinaryBigsi(bigsiArray, queryMask, numBins, subrate)

        return filteredBigsiHits
    } else {
        console.log('Query must be of length between 5Kbp and 300Kbp')
    }
}


module.exports = {
    getBloomFilterSetBitsIndices: getBloomFilterSetBitsIndices,
    getBinaryBigsiSubmatrix: getBinaryBigsiSubmatrix,
    computeQueryContainmentScores: computeQueryContainmentScores, 
    queryBinaryBigsi: queryBinaryBigsi,
    main: main
}
