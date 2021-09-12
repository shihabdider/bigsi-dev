const commonFunc = require('../common_func.js')
const { IndexedFasta, BgzipIndexedFasta } = require('@gmod/indexedfasta')

// Init a seq object (hg38 chr1) for use with tests
//const fastaPath = 'tests/data/testSeq.fa'
//const faiPath = 'tests/data/testSeq.fa.fai'
//const seq = new IndexedFasta({
//    path: fastaPath,
//    faiPath: faiPath,
//    chunkSizeLimit: 50000000
//});

test.only('makeBigsiRowPath -- should output path "test/000/001.bin"', () => {
    const result = commonFunc.makeBigsiRowPath(1, 'test')
    expect(result).toBe('test/000/001.bin')
})

test('reverseComplement - should return the correct reverse complement sequence', () => {
    const result = commonFunc.reverseComplement('ACTG')
    expect(result).toBe('CAGT')
})

test('splitIntoBuckets - should return an array with two strings corresponding to bucket regions', async () => {
    const result = await commonFunc.splitIntoBuckets(seq, '1', numBuckets=2)
    expect(result).toHaveLength(2)
    expect(result[0]).toBeDefined()
}) 

test('extractMinimizers - should return an array of hashes', async () => {
    const seqStr = await seq.getSequence('1')
    const result = commonFunc.extractMinimizers(seqStr)
    expect(result.length).toBeGreaterThan(1)
    expect(result[0]).toBeDefined()
})
