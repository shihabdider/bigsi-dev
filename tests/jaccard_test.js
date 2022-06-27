const utils = require('../bin/utils.js')
const murmur = require('murmurhash-js')
const { BloomFilter } = require('bloom-filters')

function extractKmers(sequence) {
    const kmerSize = 16
    const kmers = []
    for (let i = 0; i < (sequence.length - kmerSize + 1); i++) {
        const kmer = sequence.slice(i, i + kmerSize)
        kmers.push(kmer)
    }

    return kmers
}

function insertElementsIntoHashTable(elements) {
    const hashtable = {}
    for (const element of elements) {
        let elementHash = murmur.murmur3(element, seed=42)
        const elementNotInHashtable = !(elementHash in hashtable)
        if (elementNotInHashtable) {
            hashtable[elementHash] = true
        }
    }

    return hashtable
}

function insertElementsIntoBloomFilter(elements, bloomFilterSize){
    const numHashes = 1
    const bf = new BloomFilter(bloomFilterSize, numHashes)
    for (let element of elements){
        if (typeof element != 'string') {
            element = element.toString()
        }
        bf.add(element)
    }

    return bf
}

function getNumMatchingKmers(targetHashtable, queryKmers) {
    let numMatching = 0
    for (let element of queryKmers) {
        if (typeof element != 'string') {
            element = element.toString()
        }
        const elementHash = murmur.murmur3(element, seed=42)
        if (elementHash in targetHashtable) {
            numMatching += 1
        }
    }

    return numMatching
}

function getNumMatchingMinimizers(targetBloomFilter, queryElements) {
    let numMatching = 0
    for (let element of queryElements) {
        if (typeof element != 'string') {
            element = element.toString()
        }
        if (targetBloomFilter.has(element)) {
            numMatching += 1
        }
    }

    return numMatching
}

function computeJaccardContainment(elementsShared, elementsInQuery) {
    return elementsShared/elementsInQuery
}

function computeBloomFilterSize(n, p) {

    return Math.ceil((n * Math.log(p)) / Math.log(1 / Math.pow(2, Math.log(2))));
}

function computeJaccardDiff(targetHashtable, targetWinnowedBloomFilter, querySeq, windowSize) {

    const queryKmers = extractKmers(querySeq)
    const numMatching = getNumMatchingKmers(targetHashtable, queryKmers)
    const numInQuery = queryKmers.length 
    const true_j = computeJaccardContainment(numMatching, numInQuery)

    const queryMinimizers = utils.extractMinimizers(querySeq, windowSize)
    const numMinimizersMatching = getNumMatchingMinimizers(targetWinnowedBloomFilter, queryMinimizers)
    const numMinimizersInQuery = queryMinimizers.length
    const winnowed_j = computeJaccardContainment(numMinimizersMatching, numMinimizersInQuery)
    const jaccard_diff = true_j - winnowed_j

    return jaccard_diff
}


async function main() {
    const { argv } = require('yargs')
        .scriptName('jaccard-test')
        .usage('Usage: $0 --ref --query --window')
        .option('ref', {
            alias: 'r',
            describe: 'Fasta file of target sequence',
            demandOption: 'Fasta file is required',
            type: 'string',
            nargs: 1,
        })
        .option('query', {
            alias: 'q',
            describe: 'Fasta file of query sequence',
            demandOption: 'Fasta file is required',
            type: 'string',
            nargs: 1,
        })
        .option('windowSize', {
            alias: 'w',
            describe: 'Window size for winnowing',
            demandOption: 'Window size is required',
            type: 'number',
            nargs: 1,
        })

    const targetFai = `${argv.ref}.fai`
    const target = await utils.loadFasta(argv.ref, targetFai)
    const targetSeqNames = await target.getSequenceList()
    const targetSeq = await target.getSequence(targetSeqNames[0])

    const targetKmers = extractKmers(targetSeq)
    const targetHashtable = insertElementsIntoHashTable(targetKmers)

    const targetMinimizers = utils.extractMinimizers(targetSeq, argv.windowSize)
    const bloomFilterSizeMinimizers = computeBloomFilterSize(targetMinimizers.length, p = 0.1749)
    const targetWinnowedBloomFilter = insertElementsIntoBloomFilter(targetMinimizers, bloomFilterSizeMinimizers)
    //const targetWinnowedHashtable = insertElementsIntoHashTable(targetMinimizers)

    const queryFai = `${argv.query}.fai`
    const query = await utils.loadFasta(argv.query, queryFai)
    const querySeqNames = await query.getSequenceList()
    const querySeqs = []
    for (const seqName of querySeqNames) {
        const querySeq = await query.getSequence(seqName)
        querySeqs.push(querySeq)
    }

    const jaccard_diffs = []
    for (const querySeq of querySeqs) {
        const jaccard_diff = computeJaccardDiff(targetHashtable, targetWinnowedBloomFilter, querySeq, argv.windowSize)
        jaccard_diffs.push(jaccard_diff)
        console.log(jaccard_diff)
    }
}

main()
