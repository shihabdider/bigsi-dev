const murmur = require('murmurhash-js')
const { BloomFilter } = require('bloom-filters')
const { IndexedFasta } = require('@gmod/indexedfasta')
const fs = require('fs')

function zeroPadBitstring(bitstring, places){
    const paddedString = bitstring.padStart(places, '0')
    return paddedString
}

function reverseComplement(sequence){
    var reverseSeq=sequence.split('').reverse().join('')

    let COMPLEMENT_BASES = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G'},
        re = /[ATCG]/g;

    var revComplementSeq = reverseSeq.replace(re, function (sequence) {
        return COMPLEMENT_BASES[sequence]
    });

    return revComplementSeq
}

async function loadFasta(fastaPath, faiPath){
    const seq = new IndexedFasta({
        path: fastaPath,
        faiPath: faiPath,
        chunkSizeLimit: 5e8
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

// Based on MashMap's winnowing algorithm
function extractMinimizers(seq, windowSize){
    seq = seq.toUpperCase()

    const kmerSize = 16
    const seed = 42

    let minimizers = []
    let deque = [] // array of {hash, offset}
    let revSequence = reverseComplement(seq)
    for (let i = 0; i < (seq.length - kmerSize + 1); i++){
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

function makeMinimizersBloomFilter(minimizers, bloomFilterSize) {
    // adjust filter size based on number of inserted elements and desired false pos 
    // rate
    const numHashes = 1
    const minimizersBloomFilter = new BloomFilter(bloomFilterSize, numHashes)
    for (const minimizer of minimizers){
        minimizersBloomFilter.add(minimizer.toString())
    }

    const minimizersBloomFilterArray = minimizersBloomFilter._filter

    return minimizersBloomFilterArray
}

function writeToJSON(object, filename){
    const json = JSON.stringify(object, null, 4);

    fs.writeFile(filename, json, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`File written to: ${filename}`);
        }
    });
}

module.exports = {
    zeroPadBitstring: zeroPadBitstring,
    reverseComplement: reverseComplement,
    loadFasta: loadFasta,
    getFilteredGenomeSeqs: getFilteredGenomeSeqs,
    extractMinimizers: extractMinimizers,
    makeMinimizersBloomFilter: makeMinimizersBloomFilter,
    writeToJSON: writeToJSON,
}

