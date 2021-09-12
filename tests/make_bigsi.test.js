const makeBigsi = require('../make_bigsi.js')
const { IndexedFasta } = require('@gmod/indexedfasta')

// Init a seq object (hg38 chr1 ~30Kb subseq) for use with tests
const fastaPath = 'tests/data/testRefSeq.fa'
const faiPath = 'tests/data/testQuerySeq.fa.fai'

const SEQ = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});


// Init a seq object (hg38 3GB) for use with tests
const genomeFastaPath = 'tests/data/GCF_000001405.26_GRCh38_genomic.fna'
const genomeFaiPath = 'tests/data/GCF_000001405.26_GRCh38_genomic.fna.fai'

const GENOME = new IndexedFasta({
    path: genomeFastaPath,
    faiPath: genomeFaiPath,
    chunkSizeLimit: 50000000
});

test('makeBucketMap -- should output { object } with ref name, start and end', async () =>{
    const seqSizeThreshold = 1*10**7
    const numBuckets = 10

    const result = await makeBigsi.makeBucketToPosition(GENOME, seqSizeThreshold, numBuckets)
    expect(typeof result['0']['refName']).toBe('string')
})

test('getFilteredGenomeSeqs -- should output an array of {string} sequence names', async () =>{
    const result = await makeBigsi.getFilteredGenomeSeqs(GENOME)
    expect(typeof result[0].slice(1,2)).toBe('string')
})

test('splitSeqIntoBuckets -- should output an array of {string} sequences', async () => {
    const result = await makeBigsi.splitSeqIntoBuckets(SEQ, seqName='1', numBuckets=2)
    expect(typeof result[0].slice(1,2)).toBe('string')
})

test('makeBucketsBloomFilters -- should output an array of {int} bloom filter', async () => {
    const bucketSequences = await makeBigsi.splitSeqIntoBuckets(SEQ, seqName='1', numBuckets=2)
    const result = makeBigsi.makeBucketsBloomFilters(bucketSequences)
    expect(typeof result[0]._filter[0]).toBe('number')
})


test('buildBigsi -- should output an array of arrays of {int} bigsi matrix', async () => {
    const bucketSequences = await makeBigsi.splitSeqIntoBuckets(SEQ, seqName='1', numBuckets=2)
    const bigsiMatrix = await makeBigsi.buildBigsi(bucketSequences)
    console.log(bigsiMatrix)
    const result = bigsiMatrix()
    expect(typeof result[0][0]).toBe('number')
})

const TIMEOUT = 3*60*60*1000    //max timeout for async (3 hrs)
test('makeGenomeBigsis -- should output array of bigsi matrices for chr 1 and 2', async () => {
    const seqSizeThreshold = 2*10**8    // to get only chr 1 and 2 for test
    const result  = await makeBigsi.makeGenomeBigsis(GENOME, seqSizeThreshold)
    expect(typeof result[0]).toBe('function')
}, TIMEOUT)

test('mergeBigsis -- should output merged bigsi matrix', async () => {
    const seqSizeThreshold = 2*10**8    // to get only chr 1 and 2 for test
    const bigsis  = await makeBigsi.makeGenomeBigsis(GENOME, seqSizeThreshold)
    const result  = await makeBigsi.mergeBigsis(bigsis)
    expect(result.size()[1]).toBe(20)
}, TIMEOUT)

test('makeHexBigsi -- should output an array of {string} hexstrings for each row of bigsi', async () => {
    const seqSizeThreshold = 2*10**8    // to get only chr 1 and 2 for test
    const bigsis  = await makeBigsi.makeGenomeBigsis(GENOME, seqSizeThreshold)
    const bigsi  = await makeBigsi.mergeBigsis(bigsis)
    const result = makeBigsi.makeHexBigsi(bigsi)
    expect(typeof result[0]).toBe('string')
}, TIMEOUT)

test('makeBitBigsi -- should output array of binary { strings } for each row of bigsi', async () => {
    const seqSizeThreshold = 2*10**8    // to get only chr 1 and 2 for test
    const bigsis  = await makeBigsi.makeGenomeBigsis(GENOME, seqSizeThreshold)
    const bigsi  = await makeBigsi.mergeBigsis(bigsis)
    const result = makeBigsi.makeBitBigsi(bigsi)
    expect(typeof result[0]).toBe('string')
}, TIMEOUT)

test('writeBigsiArrayBuf -- should output binary file for each row of bigsi', async () => {
    const seqSizeThreshold = 2*10**8    // to get only chr 1 and 2 for test
    const bigsis  = await makeBigsi.makeGenomeBigsis(GENOME, seqSizeThreshold)
    const bigsi  = await makeBigsi.mergeBigsis(bigsis)
    const bitBigsi = makeBigsi.makeBitBigsi(bigsi)
    const numDirs = 300
    const root = 'test'
    const result = makeBigsi.writeBigsiArrayBuf(numDirs, root, bitBigsi)
}, TIMEOUT)
