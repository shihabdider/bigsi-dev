const murmur = require('murmurhash-js')
const { BloomFilter } = require('bloom-filters')
const { IndexedFasta } = require('@gmod/indexedfasta')
const cdf = require('binomial-cdf');

function zeroPadBitstring(bitstring, places){
    const paddedString = bitstring.padStart(places, '0')
    return paddedString
}

function zeroPad(num, places){
    const paddedString = String(num).padStart(places, '0')
    return paddedString
}

async function loadFasta(fastaPath, faiPath){
    const seq = new IndexedFasta({
        path: fastaPath,
        faiPath: faiPath,
        chunkSizeLimit: 5e7
    });

    return seq
}

/**
 * @param { IndexedFasta } genome -  sequence object for the genome
 * @param { number } seqSizeThreshold - minimum size of sequence to filter
 *
 * @returns { array of strings } seqNames - filtered seq names
 */
async function getFilteredGenomeSeqs(genome, seqSizeThreshold=10**7){
    const seqNames = await genome.getSequenceList()

    const filteredSeqNames = []
    for (let i=0; i < seqNames.length; i++){
        const seqSize = await genome.getSequenceSize(seqNames[i])
        if (seqSize >= seqSizeThreshold){
            filteredSeqNames.push(seqNames[i])
        }
    }

    return filteredSeqNames
    
}

function reverseComplement(sequence){
    let reverseSeq = sequence.split('').reverse().join('')

    const COMPLEMENT_BASES = {
        'A': 'T', 
        'T': 'A', 
        'G': 'C', 
        'C': 'G'
    }
    const re = /[ATCG]/g;

    var revComplementSeq = reverseSeq.replace(re, function (sequence) {
        return COMPLEMENT_BASES[sequence]
    });

    return revComplementSeq
}

// Based on MashMap's winnowing algorithm
function extractMinimizers(seq){
    seq = seq.toUpperCase()

    const kmerSize = 16
    const windowSize = 100
    const seed = 42

    let minimizers = []
    let deque = [] // array of {hash, offset}
    let revSequence = reverseComplement(seq)
    for (i = 0; i < (seq.length - kmerSize + 1); i++){
        let currentWindowIndex = i - windowSize + 1
        let kmerHashFwd = murmur.murmur3(seq.slice(i,i+kmerSize), seed)
        let kmerHashBwd = murmur.murmur3(revSequence.slice(-(i+kmerSize), -i), seed)
        let kmerHash = Math.min(kmerHashFwd, kmerHashBwd)

        while (deque.length != 0 && deque[0].offset <= i - windowSize){
            deque.shift()
        }

        while (deque.length != 0 && deque.slice(-1)[0].hash >= kmerHash)
        {
            deque.pop()
        }

        deque.push({'hash':kmerHash, 'offset':i})

        if (currentWindowIndex >= 0){
            if ( minimizers.length == 0 || minimizers.slice(-1)[0] != deque[0].hash )
            {
                minimizers.push(deque[0].hash)
            }
        }
    }

    return minimizers
}

function computeBloomFilterFalsePosRate(numElementsInserted, bloomFilterSize){
    const numHashes = 1

    const falsePos = (1 - math.exp(
        -1*numHashes*numElementsInserted/bloomFilterSize
    ))**numHashes
    return falsePos
}

function computeFalseHitProb(falsePosRate, maxQueryMinimizers, containmentScoreThresh){
    const numMatching = maxQueryMinimizers*containmentScoreThresh
    const falseHitProb = 1 - cdf(numMatching, maxQueryMinimizers, falsePosRate)
    return false_hit_prob
}

function computeBloomFilterSize(numElementsInserted, containmentScoreThresh, totalNumBuckets){
    // initialize set parameters
    const maxQueryMinimizers = 6_000  // 300Kbp max query = 6k minimizers
    const falseHitThresh = 1e-2
    // iterate over a array size range...
    for ( let bloomFilterSize = 0; bloomFilterSize <= 1e6; bloomFilterSize += 1e3 ){
        const falsePosRate = computeBloomFilterFalsePosRate(numElementsInserted, bloomFilterSize)
        const falseHitProb = computeFalseHitProb(
            falsePosRate, 
            maxQueryMinimizers, 
            containmentScoreThresh
        )

        // accounting for all buckets in bigsi
        const falseHitProbUpper = falseHitProb*totalNumBuckets
        // break if false hit rate less than threshold and return
        if ( falseHitProbUpper <= falseHitThresh ) {
            return bloomFilterSize
        }
    }
}

function makeMinimizersBloomFilter(
    minimizers, 
    containmentScoreThresh, 
    totalNumBuckets
) {
    // get number of minimizers
    const numInserted = minimizers.length
    const bloomFilterSize = computeBloomFilterSize(
        numInserted, 
        containmentScoreThresh, 
        totalNumBuckets
    )
    // adjust filter size based on number of inserted elements and desired false pos 
    // rate
    const minimizersBloomFilter = new BloomFilter(bloomFilterSize, nbHashes=1)
    for (const minimizer of minimizers){
        minimizersBloomFilter.add(minimizer.toString())
    }
    return minimizersBloomFilter
}

module.exports = {
    zeroPadBitstring: zeroPadBitstring,
    zeroPad: zeroPad,
    loadFasta: loadFasta,
    getFilteredGenomeSeqs: getFilteredGenomeSeqs,
    extractMinimizers: extractMinimizers,
    reverseComplement: reverseComplement,
    makeMinimizersBloomFilter: makeMinimizersBloomFilter,
}
