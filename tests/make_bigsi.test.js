const makeBigsi = require('../bin/make_bigsi.js')
const { IndexedFasta } = require('@gmod/indexedfasta')

// Init a seq object (hg38 chr1 ~2.9Mb subseq) for use with tests
const fastaPath = 'tests/test_data/testRefSeq.fa'
const faiPath = 'tests/test_data/testRefSeq.fa.fai'

const seq = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});


// Init a seq object (hg38 chr1 and 2 30mb each) for use with tests
const genomeFastaPath = 'tests/test_data/testGenomeSeq.fa'
const genomeFaiPath = 'tests/test_data/testGenomeSeq.fa.fai'

const genome = new IndexedFasta({
    path: genomeFastaPath,
    faiPath: genomeFaiPath,
    chunkSizeLimit: 50000000
});

describe('Bucket Map', () => {
    test('for a valid bucket map', async () =>{
        const seqSizes = await genome.getSequenceSizes()
        const seqName = Object.keys(seqSizes)[0]
        const seqSize = seqSizes[seqName]
        const seqIdx = 0
        const numBuckets = 16

        const bucketMap = await makeBigsi.makeBucketMap(seqName, seqSize, seqIdx, numBuckets)
        expect(typeof bucketMap['0']['refName']).toBe('string')
    })

    test('for a valid bucket map for whole genome', async () => {
        const seqSizeThreshold = 2e7
        const numBucketsPerSeq = 16
        const bucketMap = await makeBigsi.makeBucketToPositionMap(genome, seqSizeThreshold, numBucketsPerSeq)
        console.log(bucketMap)
        expect(typeof bucketMap['0']['refName']).toBe('string')
    })
})

describe('Preprocessing', () => {
    const seqName = '1'
    const numBuckets = 2

    test('for valid array of bucket sequences', async () => {
        const result = await makeBigsi.splitSeqIntoBuckets(seq, seqName, numBuckets)
        expect(typeof result[0].slice(1,2)).toBe('string')
    })

    test('for valid array of bloom filters', async () => {
        const bucketSequences = await makeBigsi.splitSeqIntoBuckets(seq, seqName, numBuckets)
        const bloomFilters = makeBigsi.makeBucketsBloomFilters(bucketSequences, numBuckets)
        console.log(bloomFilters)
        expect(typeof bloomFilters[0]._filter[0]).toBe('number')
    })
})


describe('Bigsi', () => {
    const seqName = '1'
    const numBuckets = 2
    const numBucketsGenome = 16 //for 60Mb seq we want 16 buckets
    const seqSizeThreshold = 2e7

    test('for valid bigsi matrix', async () => {
        const bucketSequences = await makeBigsi.splitSeqIntoBuckets(seq, seqName, numBuckets)
        const bigsiMatrix = await makeBigsi.buildBigsi(bucketSequences)
        //console.log(bigsiMatrix)
        const result = bigsiMatrix()
        expect(typeof result[0][0]).toBe('number')
    })

    const TIMEOUT = 3*60*60*1000    //max timeout for async (3 hrs)
    test('for valid array of bigsis for multiple sequences', async () => {
        const result  = await makeBigsi.makeGenomeBigsis(genome, numBucketsGenome, seqSizeThreshold)
        expect(typeof result[0]).toBe('function')
    }, TIMEOUT)

    test('for valid merged bigsi', async () => {
        const bigsis  = await makeBigsi.makeGenomeBigsis(genome, numBucketsGenome, seqSizeThreshold)
        const result  = await makeBigsi.mergeBigsis(bigsis)
        expect(result.size()[1]).toBe(numBucketsGenome)
    }, TIMEOUT)

    test('outputs bigsi as a flat array of 16bit ints', async () => {
        const bucketSequences = await makeBigsi.splitSeqIntoBuckets(seq, seqName, numBuckets)
        const bigsi = await makeBigsi.buildBigsi(bucketSequences)
        const intSize = 16
        const result = makeBigsi.bigsiToInts(bigsi, intSize)
        expect(typeof result[0]).toBe('number')
    })

    test('for valid binary bigsi from sequence', async () => {
        const bucketSequences = await makeBigsi.splitSeqIntoBuckets(seq, seqName, numBuckets)
        const bigsi = await makeBigsi.buildBigsi(bucketSequences)
        const intSize = 16
        const intRows = makeBigsi.bigsiToInts(bigsi, intSize)
        const result = makeBigsi.makeBinaryBigsi(intRows)
        expect(typeof result).toBe('object')
    })

    test('for valid hexbigsi', async () => {
        const bigsis  = await makeBigsi.makeGenomeBigsis(genome, numBucketsGenome, seqSizeThreshold)
        const bigsi  = await makeBigsi.mergeBigsis(bigsis)
        const result = makeBigsi.makeHexBigsi(bigsi)
        expect(typeof result[0]).toBe('string')
    }, TIMEOUT)

    test('if makeFastaBigsis returns an array of bigsis', () => {
        const result = makeBigsi.makeFastaBigsis(fasta, numBuckets)
        expect(typeof result[0].toBe('function'))
    })
})

