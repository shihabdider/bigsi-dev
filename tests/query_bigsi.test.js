const queryBigsi = require('../query_bigsi.js')
const matrix = require('matrix-js')
const { IndexedFasta } = require('@gmod/indexedfasta')
const { BloomFilter } = require('bloom-filters')
const fs = require('fs')

let bigsi = require('./data/testRefBigsi.json');
bigsi = matrix(bigsi)
//let genomeBigsi = require('./data/hg38_genome.json')
//genomeBigsi = matrix(genomeBigsi)

let binaryDumpBigsiPath = './tests/data/hg38_16int_bdump.bin'
const numCols = 32

let bitBigsiPath = './tests/data/testBigsiBinaryFiles'

let hexBigsi = require('./data/testRefHexBigsi.json')

// query fasta
const fastaPath = 'tests/data/testQuerySeq.fa'
const faiPath = 'tests/data/testQuerySeq.fa.fai'
const seq = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});

// ref fasta
const refFastaPath = 'tests/data/testRefSeq.fa'
const refFaiPath = 'tests/data/testRefSeq.fa.fai'
const refSeq = new IndexedFasta({
    path: refFastaPath,
    faiPath: refFaiPath,
    chunkSizeLimit: 50000000
});

// ref fasta genome (hg38)
const genomeFastaPath = 'tests/data/GCF_000001405.26_GRCh38_genomic.fna'
const genomeFaiPath = 'tests/data/GCF_000001405.26_GRCh38_genomic.fna.fai'
const genomeSeq = new IndexedFasta({
    path: genomeFastaPath,
    faiPath: genomeFaiPath,
    chunkSizeLimit: 50000000
});

/*
 * Basic test for frag vs. no frag query:
 *  Ref hg38 chr1 and query ACADM gene (~40Kb -> 8 frags and 800 minimizers)
 * Variations:
 *  1. Exact match vs 80% identity match: Require rebuilding the bigsi with new 
 *  size parameter
 *  2. 80% identity match artificial seq: Requires adding artifical random 
 *  mutations to ACADM
 *  3. 80% identity structural variant: Requires adding artifical structural 
 *  mutation to ACADM (insert/delete <8Kb segment)
 */

test('makeQueryFragsMinimizers - should return an array of arrays of ints', async () => {
    const querySeq = await seq.getSequence('1')
    const result = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    expect(typeof result[0][0]).toBe('number')
})

test('makeQueryFragsBloomFilters - should return an array of Bloom filters with at least 1 element inserted', 
    async () => {
        const querySeq = await seq.getSequence('1')
        const fragMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
        const result = queryBigsi.makeQueryFragsBloomFilters(fragMinimizers)
        expect(result[0]._filter).toContain(1)
})

test('getQueryBloomFilterOnesIndices - should return 1 value indices for a bloom filter', () => {
    const testBloomFilter = [1, 0, 1, 0, 0]
    const result = queryBigsi.getQueryBloomFilterOnesIndices(testBloomFilter)
    expect(result).toEqual([0,2])
})

test('initBigsiHits - should return an object of objects with properties: refName, start, end, hits ', async () => {
    const refSeqName = '1'
    const refSeqSize = await refSeq.getSequenceSize(refSeqName)
    const result = queryBigsi.initBigsiHits(refSeqSize, refSeqName, numBuckets=10, overhang=30000)
    expect(result['0']).toEqual(expect.objectContaining(
        {
            refName: expect.any(String),
            start: expect.any(Number),
            end: expect.any(Number),
            hits: 0
        }
    ))
})


test('initGenomeBigsiHits - should return an array of objects', async () => {
    const result = await queryBigsi.initGenomeBigsiHits(genomeSeq)

    const outputFilename = 'hg38_bigsi_hits.json'
    const resultJSON = JSON.stringify(result, null, 4);

    fs.writeFile(outputFilename, resultJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`genomeBigsiHits written to: ${outputFilename}`);
        }
    });

    expect(result[0]['0']).toEqual(expect.objectContaining(
        {
            refName: expect.any(Number),
            start: expect.any(Number),
            end: expect.any(Number),
            hits: 0
        }
    ))
})
test('getBinaryDumpBigsiSubmatrix - should return an array of {BitSets} for submatrix', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilters = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilters[0]._filter)

    const result = await queryBigsi.getBinaryDumpBigsiSubmatrix(binaryDumpBigsiPath, queryOnesIndices, numCols)
    //console.log(result)
    expect(typeof result[0]).toBe('object')
})

test('getBitBigsiSubmatrix - should return an array of binary {strings} for submatrix', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilters = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilters[0]._filter)

    const result = await queryBigsi.getBitBigsiSubmatrix(bitBigsiPath, queryOnesIndices)
    expect(typeof result[0]).toBe('object')
})

test('computeBitSubmatrixHits - should return object of objects for bucket hits ', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilters = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilters[0]._filter)
    const submatrix = await queryBigsi.getBitBigsiSubmatrix(bitBigsiPath, queryOnesIndices)

    const bigsiHits = {}
    await queryBigsi.computeBitSubmatrixHits(submatrix, bigsiHits)
    console.log(bigsiHits)
    expect(typeof bigsiHits['0']['hits']).toBe('number')
})

test.only('queryBinaryDumpBigsi - should return object containing fragment hits greater than 0', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)

    const result = await queryBigsi.queryBinaryDumpBigsi(binaryDumpBigsiPath, queryBloomFilter, numCols)
    console.log(result)
    expect(result['0']['hits']).toBeGreaterThan(0)
})

test('queryBitBigsi - should return object containing fragment hits greater than 0', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)

    const result = await queryBigsi.queryBitBigsi(bitBigsiPath, queryBloomFilter)
    console.log(result)
    expect(result['0']['hits']).toBeGreaterThan(0)
})

test('getHexBigsiSubmatrix - should return an array of strings for submatrix', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilter)

    const result = await queryBigsi.getHexBigsiSubmatrix(hexBigsi, queryOnesIndices)
    expect(typeof result[0]).toBe('string')
})

test('computeHexSubmatrixHits - should return object of objects for bucket hits ', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilter)
    const submatrix = await queryBigsi.getHexBigsiSubmatrix(hexBigsi, queryOnesIndices)

    const genomeBigsiHits = await queryBigsi.initGenomeBigsiHits(genomeSeq)

    await queryBigsi.computeHexSubmatrixHits(submatrix, genomeBigsiHits)
    expect(typeof genomeBigsiHits['0']['hits']).toBe('number')
})

test('queryHexBigsi - should return object containing fragment hits greater than 0', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)

    const genomeBigsiHits = await queryBigsi.initGenomeBigsiHits(genomeSeq)

    const result = await queryBigsi.queryHexBigsi(hexBigsi, queryBloomFilter, genomeBigsiHits)
    expect(result['0']['hits']).toBeGreaterThan(0)
})

test('getBigsiSubmatrix - should return a submatrix ', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilter)

    const result = await queryBigsi.getBigsiSubmatrix(bigsi, queryOnesIndices)
    expect(typeof result[0][0]).toBe('number')
})


test('computeSubmatrixHits - should return object of objects for bucket hits ', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
    const queryOnesIndices = await queryBigsi.getQueryBloomFilterOnesIndices(queryBloomFilter)
    const submatrix = await queryBigsi.getBigsiSubmatrix(bigsi, queryOnesIndices)

    const refSeqName = '1'
    const refSeqSize = await refSeq.getSequenceSize(refSeqName)
    const bigsiHits = queryBigsi.initBigsiHits(refSeqSize, refSeqName, numBuckets=10, overhang=30000)


    await queryBigsi.computeSubmatrixHits(submatrix, bigsiHits)
    expect(typeof bigsiHits['0']['hits']).toBe('number')
})

test('queryBigsi - should return object containing fragment hits greater than 0', async () => {
    const querySeq = await seq.getSequence('1')
    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)

    const refSeqName = '1'
    const refSeqSize = await refSeq.getSequenceSize(refSeqName)
    const bigsiHits = queryBigsi.initBigsiHits(refSeqSize, refSeqName, numBuckets=10, overhang=30000)

    const result = await queryBigsi.queryBigsi(bigsi, queryBloomFilter, bigsiHits)
    expect(result['0']['hits']).toBeGreaterThan(0)
})

//test('queryGenomeBigsis - should return array of objects containing fragment hits greater than 0', async () => {
//    const querySeq = await seq.getSequence('1')
//    const queryFragmentsMinimizers = await queryBigsi.makeQueryFragsMinimizers(querySeq)
//    const queryBloomFilter = await queryBigsi.makeQueryFragsBloomFilters(queryFragmentsMinimizers)
//
//    const genomeBigsiHits = await queryBigsi.initGenomeBigsiHits(genomeSeq)
//
//    const result = await queryBigsi.queryGenomeBigsis(genomeBigsi, queryBloomFilter, genomeBigsiHits)
//    console.log(result)
//    expect(result[0]['0']['hits']).toBeGreaterThan(0)
//})
