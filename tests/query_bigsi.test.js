const queryBigsi = require('../bin/query_bigsi.js')
const matrix = require('matrix-js')
const { IndexedFasta } = require('@gmod/indexedfasta')
const { BloomFilter } = require('bloom-filters')
const fs = require('fs')
const fetch  = require('cross-fetch')

const binaryBigsiPath = './test_data/hg38_chr1.bin'
const bucketMapPath = './test_data/hg38_chr1_bucket_map.json'
//const numCols = 16

//let hexBigsi = require('./data/testRefHexBigsi.json')

// query fasta
//const fastaPath = 'tests/test_data/testQuerySeq.fa'
//const faiPath = 'tests/test_data/testQuerySeq.fa.fai'

const fastaPath = '../../seqs/human/ACADM_hg38chr1.fasta'
const faiPath = '../../seqs/human/ACADM_hg38chr1.fasta.fai'
const query = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});

async function getBigsiArray(){

    //const bigsiPath = 'http://localhost:3001/public/hg38_16int_bdump.bin'
    const bigsiPath = 'http://localhost:3001/public/hg38_chr1.bin'
    const response = await fetch(bigsiPath)
    const bigsiBuffer = await response.arrayBuffer()
    //const binaryBigsiPath = 'tests/test_data/hg38_chr1.bin'
    //const bigsiBuffer = fs.readFileSync(binaryBigsiPath)
    console.log('num bytes in the bigsi buffer:', bigsiBuffer.byteLength)

    let bigsiArray = new Uint16Array(bigsiBuffer);
    console.log('bigsiArray size: ', bigsiArray.length)
    
    return bigsiArray
}

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

describe('Preprocess query', () =>{
    test('if winnowQueryFragments returns an array of minimizers', async () => {
        const querySeq = await query.getSequence('1')
        const result = await queryBigsi.winnowQueryFragments(querySeq)
        expect(typeof result[0][0]).toBe('number')
    })

    test('if makeFragmentsBloomFilters returns an array of Bloom filters with at least 1 element inserted', 
        async () => {
            const querySeq = await query.getSequence('1')
            const fragMinimizers = await queryBigsi.winnowQueryFragments(querySeq)
            const bloomFilterSize = 290_000
            const result = queryBigsi.makeFragmentsBloomFilters(fragMinimizers, bloomFilterSize)
            console.log(result)
            expect(result[0]._filter).toContain(1)
    })
})


describe('Query bigsi', () => {
    describe('Helpers', () => {
        test('if getBloomFilterSetBitsIndices returns 1 value indices for a bloom filter', () => {
            const testBloomFilter = [1, 0, 1, 0, 0]
            const result = queryBigsi.getBloomFilterSetBitsIndices(testBloomFilter)
            expect(result).toEqual([0,2])
        })
    })

    describe('Submatrix handling', () => {
        async function getRowFilter(){
            const bigsiArray = await getBigsiArray()
            const numCols = 16
            const bloomFilterSize = bigsiArray.length/numCols

            const querySeq = await query.getSequence('1')
            const queryFragmentsMinimizers = await queryBigsi.winnowQueryFragments(querySeq)
            const queryBloomFilters = await queryBigsi.makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)
            console.log(queryBloomFilters.length)
            const rowFilter = await queryBigsi.getBloomFilterSetBitsIndices(queryBloomFilters[0]._filter)

            return rowFilter
        }

        test('if getBinaryBigsiSubmatrix returns a valid submatrix', async () => {

            const bigsiArray = await getBigsiArray()
            const numCols = 16
            const rowFilter = await getRowFilter()

            const result = queryBigsi.getBinaryBigsiSubmatrix(bigsiArray, rowFilter, numCols)
            expect(typeof result()[0]).toBe('object')
        })

        test('if computeSubmatrixHits returns bucket hits for a single query fragment', async () => {
            const bigsiArray = await getBigsiArray()
            const numCols = 16
            const rowFilter = await getRowFilter()
            const submatrix = queryBigsi.getBinaryBigsiSubmatrix(bigsiArray, rowFilter, numCols)

            const bigsiHits = {}
            await queryBigsi.computeSubmatrixHits(submatrix, bigsiHits, numCols)
            console.log(bigsiHits)
            expect(typeof bigsiHits['1']['hits']).toBe('number')
        })

        test('if computeQueryContainmentScores returns valid scores for submatrix', async () => {
            const bigsiArray = await getBigsiArray()
            const numCols = 16
            const rowFilter = await getRowFilter()
            const submatrix = queryBigsi.getBinaryBigsiSubmatrix(bigsiArray, rowFilter, numCols)

            const bigsiHits = {}
            await queryBigsi.computeQueryContainmentScores(submatrix, bigsiHits)
            console.log(bigsiHits)
            expect(typeof bigsiHits['1']['containment']).toBe('number')
        })
    })

    describe('Query binary bigsi', () => {
        test('if queryBinaryBigsi returns object containing fragment hits', async () => {
            const bigsiArray = await getBigsiArray()
            const numCols = 32
            const bloomFilterSize = bigsiArray.length/numCols

            const querySeq = await query.getSequence('1')
            const queryFragmentsMinimizers = await queryBigsi.winnowQueryFragments(querySeq)
            const queryBloomFilters = await queryBigsi.makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)

            const result = await queryBigsi.queryBinaryBigsi(bigsiArray, queryBloomFilters, numCols)
            console.log(result)
            expect(result).toBeDefined()
        })

        test.only('if queryBinaryBigsi returns object containing non-frag counts', async () => {
            const bigsiArray = await getBigsiArray()
            const numCols = 16
            const bloomFilterSize = bigsiArray.length*16/numCols
            console.log('bfSize bin:', bloomFilterSize)

            const querySeq = await query.getSequence('1')
            const fragmentSizeZero = 0
            const queryFragmentsMinimizers = await queryBigsi.winnowQueryFragments(querySeq, fragmentSizeZero)
            const queryBloomFilters = await queryBigsi.makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)

            const result = await queryBigsi.queryBinaryBigsi(bigsiArray, queryBloomFilters, numCols)
            console.log(result)
            expect(result).toBeDefined()
        })
    })

    describe('Query hex bigsi', () => {
        test.only('if hexBigsi returns object containing non-frag counts', async () => {
            const hexBigsi = require('./test_data/hg38_chr1_hex.json')
            const bloomFilterSize = hexBigsi.length
            console.log('bfSize hex:', bloomFilterSize)

            const querySeq = await query.getSequence('1')
            const fragmentSizeZero = 0
            const queryFragmentsMinimizers = await queryBigsi.winnowQueryFragments(querySeq, fragmentSizeZero)
            console.log(queryFragmentsMinimizers)
            const queryBloomFilters = await queryBigsi.makeFragmentsBloomFilters(queryFragmentsMinimizers, bloomFilterSize)

            const result = queryBigsi.queryHexBigsi(hexBigsi, queryBloomFilters)
            console.log(result)
            expect(result).toBeDefined()
        })
    })
})




