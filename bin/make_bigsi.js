/* Makes a BIGSI from a reference sequence
 *
 * Input:
 *  sequence: (string)
 *  number of buckets: (number)
 *
 * Output:
 *  bigsi: (matrix)
 *
 */

const utils = require('./utils.js')
const matrix = require('matrix-js')
const cdf = require('binomial-cdf');
const config = require('../bigsi.config.json')
const writeBigsi = require('./write_bigsi.js')
const quantile = require( '@stdlib/stats-base-dists-binomial-quantile' );

/**
 * @param {number[]} nums
 * @param {number} k
 * @return {number[]}
 */
var topKFrequent = function (nums, k) {
    let map = new Map();
    let res = [];
    let bucket = Array.from({ length: nums.length + 1 }, () => []); // to create unique arrays

    // storing frequency of numbers in a map
    for (let n of nums) {
        map.set(n, map.has(n) ? 1 + map.get(n) : 1);
    }

    // Poppulate the bucket with numbers in frequency
    // as the index of the bucket
    for (const [key, value] of map.entries()) {
        bucket[value].push(key);
    }

    for (let i = bucket.length - 1; i >= 0; i--) {
        if (bucket[i].length > 0) {
            for (let n of bucket[i]) {
                res.push(n);
                if (res.length === k) return res;
            }
        }
    }
};

async function buildBigsi(sequence, bloomFilterSize, bucketCoords, currentBucketNum, ignoreMinimizers){
    const seqBloomFilters = [] // Bloom filters as arrays
    for (const coord of bucketCoords){
        const ithBucketSequence = sequence.slice(coord.bucketStart, coord.bucketEnd);
        console.log(coord['bucketStart'], coord['bucketEnd'], ithBucketSequence.length)
        const bucketMinimizers = utils.extractMinimizers(ithBucketSequence, config.windowSize)
        const k = Math.ceil(bucketMinimizers.length * 0.001/100)
        if (k >= 1) {
            const ignoreMinimizersBucket = topKFrequent(bucketMinimizers, k)
            ignoreMinimizers[currentBucketNum] = ignoreMinimizersBucket
        }
        const bucketBloomFilter = utils.makeMinimizersBloomFilter(
                bucketMinimizers, 
                bloomFilterSize
        )
        seqBloomFilters.push(bucketBloomFilter)
        currentBucketNum++
    }

    let bigsi = matrix(seqBloomFilters)
    bigsi = matrix(bigsi.trans()) // transpose to make bloom filters into columns of matrix

    return bigsi

}

function estimateNumMinimizers(seqLength){
    const numMinimizers = Math.ceil(seqLength/config.windowSize * 2)
    return numMinimizers
}

function computeFalseHitProb(falsePosRate, minQueryMinimizers, errorRate){
    if (errorRate < 1){
        const containmentScoreThresh = Math.exp(-1*errorRate*config.kmer)
        const numMatching = Math.ceil(minQueryMinimizers*containmentScoreThresh)
        const falseHitProb = 1 - cdf(numMatching, minQueryMinimizers, falsePosRate)
        return falseHitProb
    } else {
        return falsePosRate**minQueryMinimizers
    }
}

function computeBloomFilterFalsePosRate(numElementsInserted, bloomFilterSize, numHashes){
    const falsePos = (1 - Math.exp(
        -1*numHashes*numElementsInserted/bloomFilterSize
    ))**numHashes
    return falsePos
}

function errorToContainment(errorRate, kmerLength) {
    const containment_score = Math.exp(-1*errorRate*kmerLength)
    return containment_score
}

function containmentToError(containmentScore, kmerLength) {
    const errorRate = -1/kmerLength * Math.log(containmentScore)
    return errorRate
}

function computeLowerBoundErrorRate(errorRate, numMinimizersInQuery,
    confidenceInterval) {
    const containmentScore = errorToContainment(errorRate, config.kmer)
    let x = quantile(confidenceInterval, numMinimizersInQuery, containmentScore)

    const lowerBoundContainmentScore = Math.min(x / numMinimizersInQuery, 1);
    const lowerBoundError = containmentToError(lowerBoundContainmentScore, config.kmer)
    return lowerBoundError
}

function computeBloomFilterSize(maxNumElementsInserted, errorRate, totalNumBuckets){
    // initialize set parameters
    const minQueryMinimizers = estimateNumMinimizers(config.minQuerySize)
    const falseHitThresh = 1e-2
    const errorRateLower = computeLowerBoundErrorRate(errorRate, minQueryMinimizers,
                                         confidenceInterval=0.999995)
    const adjustedErrorRate = errorRate + (errorRate - errorRateLower)
    // iterate over a array size range...
    for ( let bloomFilterSize = 0; bloomFilterSize <= 5e7; bloomFilterSize += 1e3 ){
        const numHashes = 1
        const falsePosRate = computeBloomFilterFalsePosRate(maxNumElementsInserted, bloomFilterSize, numHashes)
        const falseHitProb = computeFalseHitProb(
            falsePosRate, 
            minQueryMinimizers, 
            adjustedErrorRate
        )

        // accounting for all buckets in bigsi
        const falseHitProbUpper = falseHitProb*totalNumBuckets
        // break if false hit rate less than threshold and return
        if ( falseHitProbUpper <= falseHitThresh ) {
            console.log('num hashes:', numHashes)
            console.log('bf false positive rate:', falsePosRate)
            return bloomFilterSize
        }
    }
}

function computeNumBuckets(seqLength) {
    const numBuckets = Math.ceil(seqLength/config.bucketSize)
    return numBuckets
}

/* estimate bloom filter size using minimizer count computed from bucket size 
 * of longest sequence + bucket overhangs
*/ 
function estimateBloomFilterSize(seqSizes, bucketSize){
    const numElementsInserted = estimateNumMinimizers(bucketSize)
    console.log('number of minimizers: ', numElementsInserted)
    const errorRate = config.errorRate
    console.log('max error rate:', errorRate)
    let totalNumBuckets = 0
    const seqSizesArr = Object.values(seqSizes)
    for (const seqSize of seqSizesArr) {
        const numBucketsInSeq = computeNumBuckets(seqSize)
        totalNumBuckets += numBucketsInSeq
    }

    const bloomFilterSize = computeBloomFilterSize(
        numElementsInserted, 
        errorRate,
        totalNumBuckets
    )
    console.log(`optimal bloom filter size: ${bloomFilterSize}`)

    return bloomFilterSize
}


function computeBucketCoords(seqLength) {
    const bucketCoords = []
    const numBuckets = computeNumBuckets(seqLength)
    for (let bucketNum=0; bucketNum < numBuckets; bucketNum++){
        const bucketStart = Math.max(bucketNum*config.bucketSize - config.bucketOverhang, 0)
        let bucketEnd = Math.min(bucketStart + config.bucketSize + 2*config.bucketOverhang, seqLength)

        if (bucketStart === 0) { // handle first bucket 
            bucketEnd -= config.bucketOverhang
        }

        const coord = { bucketStart, bucketEnd }
        bucketCoords.push(coord)
    }

    return bucketCoords
}

/**
 * @param { IndexedFasta } fasta - indexedFasta object
 *
 * @returns { string[][] } fastaBigsis - array for each seq in fasta in 
 * bitstring format (each element is a bitstring corresponding to the bigsi row 
 * for each sequence)
 * @returns { Object } ignoreMinimizers - high freq minimizers to ignore for 
 * each bucket
 */
async function makeFastaBigsis(fasta){
    const seqNames = await fasta.getSequenceList()
    console.log('seqNames: ', seqNames)

    const seqSizes = await fasta.getSequenceSizes()
    const fullBucketSize = config.bucketSize + 2*config.bucketOverhang
    const bloomFilterSize = estimateBloomFilterSize(seqSizes, fullBucketSize)

    const fastaBigsis = []
    const ignoreMinimizers = {}
    let currentBucketNum = 0
    for (const seqName of seqNames){
        const sequence = await fasta.getSequence(seqName)
        const bucketCoords = computeBucketCoords(sequence.length, config.bucketSize)
        const seqBigsi = await buildBigsi(sequence, bloomFilterSize, bucketCoords, currentBucketNum, ignoreMinimizers)
        const bigsiBitstrings = writeBigsi.bigsiToBitstrings(seqBigsi)
        //console.log('seqBigsi', seqBigsi)
        console.log(`Bigsi of ${seqName} built...`)
        fastaBigsis.push(bigsiBitstrings)
        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process used ${memoryUsed} MB`)
    }

    return [fastaBigsis, ignoreMinimizers]
}

function mergeBigsis(bigsis){
    const mergedBigsi = bigsis.reduce((a, b) => a.map((x, i) => x + b[i]));
    const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
    console.log(`Process used ${memoryUsed} MB`)

    return mergedBigsi
}


/**
 * @param { IndexedFasta } fasta - indexedFasta object
 * @returns { string[] } bigsi - array of bitstrings corresponding to each row 
 * of the full bigsi
 */
async function main(fasta) {
    const minSeqLength = 30e6
    const seqSizes = Object.values(await fasta.getSequenceSizes())
    const areFastaSeqsValidSize = Math.min(...seqSizes) > minSeqLength 

    if (areFastaSeqsValidSize && config.bucketSize > 0) {
        const [bigsis, ignoreMinimizers] = await makeFastaBigsis(fasta)
        console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
        const bigsi = await mergeBigsis(bigsis)
        const paddingSize = config.intBits - (bigsi[0].length % config.intBits)
        const numColsPadded = bigsi[0].length + paddingSize
        const bigsiDims = { 'rows': bigsi.length, 'cols': numColsPadded, 'padding': paddingSize }
        console.log(`Bigsis merged!`)
        console.log('Number of (rows, cols):', bigsiDims)

        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        return [bigsi, bigsiDims, ignoreMinimizers]

    } else { 
        if (!areFastaSeqsValidSize) { 
            console.log('All sequences must be at least 30Mbp in length.') 
        }
        if (!(config.bucketSize > 0)) { 
            console.log('Bucket size must be greater than 0.') 
        }
    }
}

module.exports = {
    main: main
}
