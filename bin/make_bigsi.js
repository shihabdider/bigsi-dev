/* Makes a bigsi from a reference sequence
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
function estimateBloomFilterSize(seqSizes){
    const seqSizesArr = Object.values(seqSizes)
    const maxSeqLength = Math.ceil(Math.max(...seqSizesArr)/config.numBuckets + config.bucketOverhang)
    console.log('maximum sequence length: ', maxSeqLength)
    const maxNumElementsInserted = estimateNumMinimizers(maxSeqLength)
    console.log('max number of minimizers: ', maxNumElementsInserted)
    const errorRate = config.errorRate
    console.log('max error rate:', errorRate)
    const totalNumBuckets = seqSizesArr.length*config.numBuckets
    console.log('total number of buckets:', totalNumBuckets)

    const bloomFilterSize = computeBloomFilterSize(
        maxNumElementsInserted, 
        errorRate,
        totalNumBuckets
    )
    console.log(`optimal bloom filter size: ${bloomFilterSize}`)

    return bloomFilterSize
}

function computeBucketCoords(seqLength, bucketSize) {
    const bucketCoords = []
    for (let i=0; i < config.numBuckets; i++){
        const bucketStart = Math.max((bucketSize - config.bucketOverhang)*i, 0)
        let bucketEnd = bucketStart + bucketSize

        // Handle last bucket
        if (i == config.numBuckets - 1){
            bucketEnd = seqLength
        }

        const coord = { bucketStart, bucketEnd }
        bucketCoords.push(coord)
    }

    return bucketCoords
}

/**
 * @param { IndexedFasta } fasta - indexedFasta object
 *
 * @returns { matrix[] } fastaBigsis - bigsi matrix for each seq in fasta in 
 * compressed format (each column is a vector of 16bit ints corresponding to 
 * a 16-column row)
 */
async function makeFastaBigsis(fasta){
    const seqNames = await fasta.getSequenceList()
    console.log('seqNames: ', seqNames)

    const seqSizes = await fasta.getSequenceSizes()
    const bloomFilterSize = estimateBloomFilterSize(seqSizes)

    const fastaBigsis = []
    for (const seqName of seqNames){
        const sequence = await fasta.getSequence(seqName)
        const bucketSize = Math.round(
            (seqSizes[seqName] + ((config.numBuckets - 1)*config.bucketOverhang))/config.numBuckets
        )
        const bucketCoords = computeBucketCoords(sequence.length, bucketSize)
        const seqBigsi = await buildBigsi(sequence, bloomFilterSize, bucketCoords)
        //console.log('seqBigsi', seqBigsi)
        const seqBigsiInts = writeBigsi.bigsiToInts(seqBigsi)
        //console.log('seqBigsiInts', seqBigsiInts)
        const seqBigsiVector = seqBigsiInts.map((elem) => [elem])
        //console.log('seqBigsiVector', seqBigsiVector)
        console.log(`Bigsi of ${seqName} built...`)
        fastaBigsis.push(matrix(seqBigsiVector))
        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process used ${memoryUsed} MB`)
    }

    return fastaBigsis
}

function mergeBigsis(bigsis){
    let mergedBigsi = bigsis[0]
    if (bigsis.length > 1){
        for (let i=1; i < bigsis.length; i++){
            mergedBigsi = matrix(mergedBigsi.merge.right(bigsis[i]()))
        }
    }
    const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
    console.log(`Process used ${memoryUsed} MB`)

    return mergedBigsi
}

async function main(fasta) {
    const minSeqLength = 30e6
    const seqSizes = Object.values(await fasta.getSequenceSizes())
    const areFastaSeqsValidSize = Math.min(...seqSizes) > minSeqLength 

    if (areFastaSeqsValidSize && config.numBuckets > 0) {
        const bigsis = await makeFastaBigsis(fasta, config.numBuckets)
        console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
        const bigsi = await mergeBigsis(bigsis)
        console.log(`Bigsis merged!`)
        console.log('Number of (rows, cols):', bigsi.size())

        const memoryUsed = process.memoryUsage().rss / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        return bigsi
    } else { 
        if (!areFastaSeqsValidSize) { 
            console.log('All sequences must be at least 30Mbp in length.') 
        }
        if (!(config.numBuckets > 0)) { 
            console.log('Number of buckets must be greater than 0.') 
        }
    }
}

module.exports = {
    main: main
}
