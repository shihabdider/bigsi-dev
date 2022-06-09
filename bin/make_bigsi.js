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

function makeBucketBloomFilter(sequence, bloomFilterSize){
    const bucketMinimizers = utils.extractMinimizers(sequence, config.windowSize)
    const bucketBloomFilter = utils.makeMinimizersBloomFilter(
            bucketMinimizers, 
            bloomFilterSize
        )

    return bucketBloomFilter
}

async function buildBigsi(sequence, bloomFilterSize, bucketCoords){
    const seqBloomFilters = [] // Bloom filters as arrays
    for (const coord of bucketCoords){
        const ithBucketSequence = sequence.slice(coord.bucketStart, coord.bucketEnd);
        console.log(coord['bucketStart'], coord['bucketEnd'], ithBucketSequence.length)
        const bucketBloomFilter = makeBucketBloomFilter(ithBucketSequence, bloomFilterSize)
        seqBloomFilters.push(bucketBloomFilter)
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

function computeNumHashes(bloomFilterSize, numElementsInserted) {
    const numHashes = Math.ceil((bloomFilterSize/numElementsInserted)*Math.log(2))
    return numHashes
}

function computeBloomFilterFalsePosRate(numElementsInserted, bloomFilterSize, numHashes){
    const falsePos = (1 - Math.exp(
        -1*numHashes*numElementsInserted/bloomFilterSize
    ))**numHashes
    return falsePos
}

function computeBloomFilterSize(maxNumElementsInserted, errorRate, totalNumBuckets){
    // initialize set parameters
    const minQueryMinimizers = 2*config.minQuerySize/config.windowSize
    const falseHitThresh = 1e-2
    // iterate over a array size range...
    for ( let bloomFilterSize = 0; bloomFilterSize <= 5e7; bloomFilterSize += 1e3 ){
        const numHashes = 1
        const falsePosRate = computeBloomFilterFalsePosRate(maxNumElementsInserted, bloomFilterSize, numHashes)
        const falseHitProb = computeFalseHitProb(
            falsePosRate, 
            minQueryMinimizers, 
            errorRate
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

/* estimate bloom filter size using minimizer count computed from bucket size 
 * of longest sequence + bucket overhang
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

function computeNumBuckets(seqLength) {
    const numBuckets = Math.ceil(seqLength/config.bucketSize)
    return numBuckets
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
 */
async function makeFastaBigsis(fasta){
    const seqNames = await fasta.getSequenceList()
    console.log('seqNames: ', seqNames)

    const seqSizes = await fasta.getSequenceSizes()
    const fullBucketSize = config.bucketSize + 2*config.bucketOverhang
    const bloomFilterSize = estimateBloomFilterSize(seqSizes, fullBucketSize)

    const fastaBigsis = []
    for (const seqName of seqNames){
        const sequence = await fasta.getSequence(seqName)
        const bucketCoords = computeBucketCoords(sequence.length, config.bucketSize)
        const seqBigsi = await buildBigsi(sequence, bloomFilterSize, bucketCoords)
        const bigsiBitstrings = writeBigsi.bigsiToBitstrings(seqBigsi)
        //console.log('seqBigsi', seqBigsi)
        console.log(`Bigsi of ${seqName} built...`)
        fastaBigsis.push(bigsiBitstrings)
        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process used ${memoryUsed} MB`)
    }

    return fastaBigsis
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
        const bigsis = await makeFastaBigsis(fasta, config.numBuckets)
        console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
        const bigsi = await mergeBigsis(bigsis)
        const bigsiDims = [bigsi.length, bigsi[0].length]
        console.log(`Bigsis merged!`)
        console.log('Number of (rows, cols):', bigsiDims)

        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        return bigsi

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
