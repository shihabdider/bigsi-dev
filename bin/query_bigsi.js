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
 *  minimizer and is set 100, meaning each minimizer "represents" 50 k-mers.
 *  
 */

const cdf = require('binomial-cdf')
const matrix = require('matrix-js')
const commonFunc = require('./common_func.js')
const BitSet = require('bitset')
const fs = require('fs')
//import Uint1Array from 'uint1array'
const Uint1Array = require('uint1array').default;


// One function is used for fragmenting and winnowing to prevent double 
// iteration
async function winnowQueryFragments(querySeq, fragmentSize=2500){

    const querySize = querySeq.length

    const queryFragmentsMinimizers = []

    if (querySize <= 300_000 && querySize >= fragmentSize){

        for (let i=0; i < Math.floor(querySize/fragmentSize); i++){
            const start = i*fragmentSize
            const end = Math.min(start + fragmentSize, querySize)
            const fragmentSeq = querySeq.slice(start, end)
            const fragmentMinimizers = commonFunc.extractMinimizers(queryFragmentSeq)
            queryFragmentsMinimizers.push(fragmentMinimizers)
        }
    } else { console.error("Inappropriate query size") }

    return queryFragmentsMinimizers
}

function makeFragmentsBloomFilters(queryFragmentsMinimizers){

    const fragmentsBloomFilters = []
    for (let fragmentMinimizerSet of queryFragmentsMinimizers){
        const bf = commonFunc.makeMinimizersBloomFilter(fragmentMinimizerSet)
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

function initBigsiHits(seqSize, seqName, numBuckets=16, bucketOverhang=300_000){
    const bucketSize = Math.round(seqSize/(numBuckets-1))

    const bigsiHits = {} 
    for (let bucketNum=0; bucketNum < numBuckets; bucketNum++){
        bigsiHits[bucketNum] = {}

        if (bucketNum == 0){
            bigsiHits[bucketNum] = Object.assign(
                {
                    refName: seqName,
                    start: 0,
                    end: bucketSize,
                    hits: 0,
                    score: 0,
                }, 
                bigsiHits[bucketNum])
        } else {
            const intervalStart = (bucketNum*(bucketSize-bucketOverhang))+1
            let intervalEnd = intervalStart + bucketSize

            if (bucketNum == numBuckets - 1){
                intervalEnd = seqSize
            }

            bigsiHits[bucketNum] = Object.assign(
                {
                    refName: seqName,
                    start: intervalStart,
                    end: intervalEnd,
                    hits: 0,
                    score: 0,
                }, 
                bigsiHits[bucketNum]
            )
        }
    }

    return bigsiHits
}

/**
 * @param { IndexedFasta} genome - IndexedFasta seq object of genome
 *
 * @returns { array of objects } genomeBigsiHits - array of empty objects for storing query hits
 */
async function initGenomeBigsiHits(genome){
    const filteredSeqNames = await commonFunc.getFilteredGenomeSeqs(genome)
    const genomeBigsiHits = []
    for (let i=0; i < filteredSeqNames.length; i++){
        const seqName = filteredSeqNames[i]
        const seqSize = await genome.getSequenceSize(seqName)
        const seqBigsiHits = initBigsiHits(seqSize, parseInt(i+1))

        genomeBigsiHits.push(seqBigsiHits)
    }

    return genomeBigsiHits
}


/** Converts a buffer to a u16 typed array
 */
function toArrayBuffer(buffer) {
    const buf = new ArrayBuffer(buffer.length);
    const view = new Uint16Array(buf);
    for (let i = 0; i < buffer.length; ++i) {
        view[i] = buffer[i];
    }
    return buf;
}

/**
 * @param { string } filename - name of binary dump bigsi file
 * @param { array of numbers } rowFilter - array of row numbers for query rows
 * @param { number } numCols - number of columns in bigsi
 *
 * @returns { array of strings } submatrix - array of bitsets of query rows 
 * submatrix
 */
async function getBinaryDumpBigsiSubmatrix(filename, rowFilter, numCols){
    const submatrix = []

    try {
        let bigsiRowBuf = await fs.readFileSync(filename, (err, data) => {
                if (err) {
                    console.error(err)
                    return
                }
            })

        const array = new Uint16Array(toArrayBuffer(bigsiRowBuf));
        console.log('array size: ', array.length)
        const numSeqs = numCols/16

        for (let rowNum of rowFilter){

            const offsetStart = rowNum*numSeqs
            const offsetEnd = offsetStart + numSeqs
            console.log('offsets: ', offsetStart, offsetEnd)
            console.log('ints: ', array.slice(offsetStart, offsetEnd))

            const rowBitString = Array.from(array.slice(offsetStart, offsetEnd))
                .map((num) => {return num.toString(2)})
                .join('')
            console.log('bitstring: ', rowBitString)
            const bs = new BitSet(rowBitString)
            submatrix.push(bs);
        }
    } catch(err) {
        console.log(err)
    }


    return submatrix
}

/**
 * @param { string } root - parent directory where bigsi is stored
 * @param { array of numbers } rowFilter - array of row numbers for query rows
 *
 * @returns { array of strings } submatrix - array of bitsets of query rows 
 * submatrix
 */
async function getBitBigsiSubmatrix(root, rowFilter){
    const submatrix = []

    for (let i = 0; i < rowFilter.length; i++){
        const rowPath = commonFunc.makeBigsiRowPath(rowFilter[i], root)

        try {
            let bigsiRowBuf = await fs.promises.readFile(rowPath, (err, data) => {
                    if (err) {
                        console.error(err)
                        return
                    }
                })

            const array = new Uint8Array(bigsiRowBuf);
            const bs = new BitSet(array)
            submatrix.push(bs);

        } catch(err) {
            console.log(err)
        }

    }

    return submatrix
}

/**
 * @param { array of strings } submatrix - array of bitsets for query rows 
 * of bigsi
 * @param { bigsiHits } bigsiHits - empty object for storing returned hits 
 *
 * @returns - bigsiHits with updated hits attributes 
 */
function computeBitSubmatrixHits(submatrix, bigsiHits){
    let bucketHits = (new BitSet).flip() // init bitset to all 1s for bitwise AND
    for (let i = 0; i < submatrix.length; i++){
        bucketHits = bucketHits.and(submatrix[i])
    }

    const hitsBucketNums = bucketHits.toArray().map( value => bucketHits.toString().length - 1 - value )

    for (let i = 0; i < hitsBucketNums.length; i++){
        const bucketNum = hitsBucketNums[i]
        if (bucketNum in bigsiHits){
            bigsiHits[bucketNum]['hits'] += 1
        } else {
            bigsiHits[bucketNum] = {'hits': 1}
        }

    }
}

/** 
 * @param { string } filename - name of binary dump bigsi file
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { string } numCols - number of columns in bigsi
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryBinaryDumpBigsi(filename, queryFragmentsBloomFilters, numCols){

    const bigsiHits = {}

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query BFs: ', queryFragmentsBloomFilters.length)

    for (let m = 0; m < queryFragmentsBloomFilters.length; m++){
        const queryBFOnesIndices = getBloomFilterSetBitsIndices(queryFragmentsBloomFilters[m]._filter)
        
        let querySubmatrix = await getBinaryDumpBigsiSubmatrix(filename, queryBFOnesIndices, numCols)
        computeBitSubmatrixHits(querySubmatrix, bigsiHits)
    }

    for (let bucketId in bigsiHits) {
        bigsiHits[bucketId]['score'] = `${bigsiHits[bucketId]['hits']}/${numFragments}`;
    }

    return bigsiHits
}

/** 
 * @param { string } root - parent directory where bigsi rows are stored
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { object } bigsiHits - empty object for storing query hits
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryBitBigsi(root, queryFragmentsBloomFilters){

    const bigsiHits = {}

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query BFs: ', queryFragmentsBloomFilters.length)

    for (let m = 0; m < queryFragmentsBloomFilters.length; m++){
        const queryBFOnesIndices = getBloomFilterSetBitsIndices(queryFragmentsBloomFilters[m]._filter)
        
        let querySubmatrix = await getBitBigsiSubmatrix(root, queryBFOnesIndices)
        computeBitSubmatrixHits(querySubmatrix, bigsiHits)
    }

    for (let bucketId in bigsiHits) {
        bigsiHits[bucketId]['score'] = `${bigsiHits[bucketId]['hits']}/${numFragments}`;
    }

    return bigsiHits
}

/**
 * @param { array of strings } hexBigsi - array of hexstrings for each row of 
 * bigsi
 *
 * @returns { array of strings } submatrix - array of hexstrings of query rows 
 * submatrix
 */
function getHexBigsiSubmatrix(hexBigsi, rowFilter){
    let submatrix = []

    rowFilter.forEach(i => submatrix.push(hexBigsi[i]));

    return submatrix

}

/**
 * @param { array of strings } submatrix - array of hexstrings for query rows 
 * of bigsi
 * @param { bigsiHits } bigsiHits - empty object for storing returned hits 
 */
function computeHexSubmatrixHits(submatrix, bigsiHits){
    let bucketHits = (new BitSet).flip() // init bitset to all 1s for bitwise AND
    for (let i = 0; i < submatrix.length; i++){
        const bs = BitSet.fromHexString(submatrix[i])
        bucketHits = bucketHits.and(bs)
    }

    const bucketIds = bucketHits.toArray()
    console.log(bucketIds)
    console.log(bigsiHits)
    bucketIds.forEach(bucket => bigsiHits[bucket]['hits'] += 1)

}
    
/** 
 * @param { array of strings } bigsi - bigsi matrix as array of hexstrings
 * @param { array of bloom filters } queryFragmentsBloomFilters - an array of Bloom filters with fragment 
 * minimizers inserted
 * @param { object } bigsiHits - empty object for storing query hits
 *
 * @return { object } filteredBigsiHits - object containing fragment hits in the bigsi buckets
 */
async function queryHexBigsi(bigsi, queryFragmentsBloomFilters, bigsiHits){

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query BFs: ', queryFragmentsBloomFilters.length)

    for (let m = 0; m < queryFragmentsBloomFilters.length; m++){
        const queryBFOnesIndices = getBloomFilterSetBitsIndices(queryFragmentsBloomFilters[m]._filter)
        
        const querySubmatrix = getHexBigsiSubmatrix(bigsi, queryBFOnesIndices)
        computeHexSubmatrixHits(querySubmatrix, bigsiHits)
    }

    for (let bucketId in bigsiHits) {
        if (bigsiHits.hasOwnProperty(bucketId)) {
            bigsiHits[bucketId]['score'] = bigsiHits[bucketId]['hits']/numFragments;
        }
    }

    const filteredBigsiHits = Object.fromEntries(
        Object.entries(bigsiHits).filter(([key, value]) => value.score > 0) 
    )

    return filteredBigsiHits
}

function getBigsiSubmatrix(bigsi, rowFilter){
    /* Inputs:
     *  bigsi -- bigsi matrix object
     *  rowFilter -- array of row indices to be extracted
     *
     * Outputs:
     *  submatrix -- resulting submatrix from extracting rows
     *
     */
    let submatrix = []
    for (i = 0; i < rowFilter.length; i++){
        const row = bigsi(rowFilter[i])
        submatrix.push(row)
    }

    submatrix = matrix(submatrix).trans()

    return submatrix

}

function computeSubmatrixHits(submatrix, bigsiHits, BFFalsePositiveProb=0.82, hitFalsePositiveProb=0.001){
    /* Inputs:
     *  submatrix -- bigsi array of arrays of hit rows
     *  bigsiHits -- object to store hits results for each bucket in bigsi
     *  BFFalsePositiveProb -- probability of a false positive for a Bloom Filter 
     *      query
     *  hitFalsePositiveProb -- set probability for the false positive of a hit
     *
     * Outputs:
     *  bigsiHits -- with updated hits property for each bucket
     *
     */
    const thresholdProbability = 1 - hitFalsePositiveProb

    for (let bucketNum = 0; bucketNum < submatrix.length; bucketNum++){
        rowsum = submatrix[bucketNum].reduce((a, b) => a + b, 0)
        const score = cdf(rowsum, submatrix[bucketNum].length, BFFalsePositiveProb)

        if (score >= thresholdProbability){
            bigsiHits[bucketNum]['hits'] += 1
        }
    }

}

async function queryBigsi(bigsi, queryFragmentsBloomFilters, bigsiHits){
    /* Inputs:
     *  bigsi -- bigsi matrix
     *  queryFragmentsBloomFilters -- an array of Bloom filters with fragment 
     *  minimizers inserted
     *
     * Outputs:
     *  filteredBigsiHits -- object containing fragment hits in the bigsi buckets
     *  
     */

    const numFragments = queryFragmentsBloomFilters.length
    console.log('number of query BFs: ', queryFragmentsBloomFilters.length)

    for (let m = 0; m < queryFragmentsBloomFilters.length; m++){
        const queryBFOnesIndices = getBloomFilterSetBitsIndices(queryFragmentsBloomFilters[m]._filter)
        
        const querySubmatrix = getBigsiSubmatrix(bigsi, queryBFOnesIndices)
        computeSubmatrixHits(querySubmatrix, bigsiHits)
    }

    for (let bucketId in bigsiHits) {
        if (bigsiHits.hasOwnProperty(bucketId)) {
            bigsiHits[bucketId]['score'] = bigsiHits[bucketId]['hits']/numFragments;
        }
    }

    const filteredBigsiHits = Object.fromEntries(
        Object.entries(bigsiHits).filter(([key, value]) => value.score > 0) 
    )

    return filteredBigsiHits
}

/**
 * @param { function } genomeBigsi - matrix containing bigsi of entire genome
 * @param { array of objects } genomeBigsiHits - array of empty objects for 
 *  storing query hits
 * @param { number } numBuckets - number of buckets
 *
 * @returns { array of objects } filteredGenomeBigsiHits - array of filtered 
 *  objects containing query hits 
 */
async function queryGenomeBigsis(genomeBigsi, queryFragmentsBloomFilters, genomeBigsiHits, numBuckets=10){
    let filteredGenomeBigsiHits = []
    for (let i=0; i < genomeBigsiHits.length; i++){
        const bigsiStartColumn = i*numBuckets
        const bigsiEndColumn = (bigsiStartColumn + numBuckets - 1)
        const seqBigsi = matrix(genomeBigsi([],[bigsiStartColumn,bigsiEndColumn]))
        const seqBigsiHits = genomeBigsiHits[i]
        const filteredSeqHits = await queryBigsi(seqBigsi, queryFragmentsBloomFilters, seqBigsiHits, numBuckets)
        filteredGenomeBigsiHits.push(filteredSeqHits)
    }

    return filteredGenomeBigsiHits
}

// main is only called for tests
async function main(bigsi, querySeq){
    const queryFragmentsMinimizers = await winnowQueryFragments(querySeq)
    const queryBloomFilter = await makeFragmentsBloomFilters(queryFragmentsMinimizers)

    const refSeqName = '1'
    const refSeqSize = await refSeq.getSequenceSize(refSeqName)
    const bigsiHits = initBigsiHits(refSeqSize, refSeqName, numBuckets=10, bucketOverhang=30000)

    //const filteredBigsiHits = await queryBigsi(bigsi, queryBloomFilter, bigsiHits)
    //const filteredBigsiHits = await queryHexBigsi(bigsi, queryBloomFilter, bigsiHits)
    
    const root = 'bigsiBinaryFiles'
    const filteredBigsiHits = await queryBitBigsi(root, queryBloomFilter, bigsiHits)

    return filteredBigsiHits
}

module.exports = {
    winnowQueryFragments: winnowQueryFragments,
    makeFragmentsBloomFilters: makeFragmentsBloomFilters,
    getBloomFilterSetBitsIndices: getBloomFilterSetBitsIndices,
    initBigsiHits: initBigsiHits,
    initGenomeBigsiHits: initGenomeBigsiHits,
    getBinaryDumpBigsiSubmatrix: getBinaryDumpBigsiSubmatrix,
    getBitBigsiSubmatrix: getBitBigsiSubmatrix,
    computeBitSubmatrixHits: computeBitSubmatrixHits,
    queryBinaryDumpBigsi: queryBinaryDumpBigsi,
    queryBitBigsi: queryBitBigsi,
    getHexBigsiSubmatrix: getHexBigsiSubmatrix,
    computeHexSubmatrixHits: computeHexSubmatrixHits,
    queryHexBigsi: queryHexBigsi,
    getBigsiSubmatrix: getBigsiSubmatrix,
    computeSubmatrixHits: computeSubmatrixHits,
    queryBigsi: queryBigsi,
    queryGenomeBigsis: queryGenomeBigsis
}
