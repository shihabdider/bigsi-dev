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

function makeSeqBloomFilter(sequence, bloomFilterSize){
    const seqMinimizers = utils.extractMinimizers(sequence, config.windowSize)
    const seqBloomFilter = utils.makeMinimizersBloomFilter(
            seqMinimizers, 
            bloomFilterSize
        )

    return seqBloomFilter
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
    console.log('adjusted error rate: ', adjustedErrorRate)
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
            return bloomFilterSize
        }
    }
}

function computeNumBuckets(seqLength) {
    const numBuckets = Math.ceil(seqLength/config.bucketSize)
    return numBuckets
}

/* estimate bloom filter size using minimizer count computed from largest 
 * sequencee
*/ 
function estimateBloomFilterSize(seqSizes){
    const seqSizesArr = Object.values(seqSizes)
    const numElementsInserted = estimateNumMinimizers(Math.max(...seqSizesArr))
    console.log('number of minimizers: ', numElementsInserted)
    const errorRate = config.errorRate
    console.log('max error rate:', errorRate)

    const totalNumBuckets = seqSizesArr.length
    const bloomFilterSize = computeBloomFilterSize(
        numElementsInserted, 
        errorRate,
        totalNumBuckets
    )
    console.log(`optimal bloom filter size: ${bloomFilterSize}`)

    return bloomFilterSize
}

/**
 * @param { IndexedFasta } fasta - indexedFasta object
 *
 * @returns { string[] } fastaBigsi - array for each seq in fasta in 
 * bitstring format (each element is a bitstring corresponding to the bigsi row 
 * for each sequence)
 */
async function makeFastaBigsi(fasta){
    const seqNames = await fasta.getSequenceList()
    console.log('seqNames: ', seqNames)

    const seqSizes = await fasta.getSequenceSizes()
    const bloomFilterSize = estimateBloomFilterSize(seqSizes)

    const seqBloomFilters = []
    for (const seqName of seqNames){
        const sequence = await fasta.getSequence(seqName)
        const seqBloomFilter = makeSeqBloomFilter(sequence, bloomFilterSize)
        seqBloomFilters.push(seqBloomFilter)
        //console.log('seqBigsi', seqBigsi)
        console.log(`Bigsi of ${seqName} built...`)
        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process used ${memoryUsed} MB`)
    }

    let bigsi = matrix(seqBloomFilters)
    bigsi = matrix(bigsi.trans()) // transpose to make bloom filters into columns of matrix

    const fastaBigsi = writeBigsi.bigsiToBitstrings(bigsi)
    return fastaBigsi
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

    if (areFastaSeqsValidSize) {
        const bigsi = await makeFastaBigsi(fasta)
        console.log(`Bigsi for ${bigsi.length} sequences created, merging...`)
        const bigsiDims = { 'rows': bigsi.length, 'cols': seqSizes.length }
        console.log('Number of (rows, cols):', bigsiDims)

        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        return [bigsi, bigsiDims]

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
