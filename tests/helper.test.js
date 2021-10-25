const helper = require('../bin/helper')
const { IndexedFasta, BgzipIndexedFasta } = require( '@gmod/indexedfasta')

// Init a seq object (hg38 chr1) for use with tests
const fastaPath = 'tests/test_data/testRefSeq.fa'
const faiPath = 'tests/test_data/testRefSeq.fa.fai'
const testSeq = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});

describe('General helper functions', () => {
    it('pads the bitstring to 8 bits', () => {
        const bitstring = '01010'
        const result = helper.zeroPadBitstring(bitstring, 8)
        expect(result).toBe('00001010')
    })


    it('returns the correct reverse complement sequence', () => {
        const result = helper.reverseComplement('ACTG')
        expect(result).toBe('CAGT')
    })
})

describe('Sequence Handling', () => {
    it('loads a fasta file into a IndexedFasta', () => {
        // to-do
    })

    it('prints a list of names of genome sequences filtered by a size threshold', () => {
        // to-do
    })
})

describe('Bloom Filter', () => {
    it('returns an array of hashes of minimizers', async () => {
        const seq = await testSeq.getSequence('1')
        const result = helper.extractMinimizers(seq)
        expect(result.length).toBeGreaterThan(1)
        expect(result[0]).toBeDefined()
    })

    it('returns false positive rate for a Bloom filter with given parameters', () => {
        const numElementsInserted = 300_000
        const bloomFilterSize = 312_000
        const falsePosRate = helper.computeBloomFilterFalsePosRate(numElementsInserted, bloomFilterSize)
        console.log(falsePosRate)
        expect(falsePosRate).toBeCloseTo(0.618)
    })

    it('returns false hit probability', () => {
        const falsePosRate = 0.618
        const minQueryMinimizers = 100
        const containmentScoreThresh = 0.8
        const falseHitProb = helper.computeFalseHitProb(falsePosRate, minQueryMinimizers, containmentScoreThresh)
        expect(falseHitProb).toBeCloseTo(2.723e-05)
    })

    it('returns optimal size of a Bloom filter with given parameters', () => {
        const numElementsInserted = 300_000
        const containmentScoreThresh = 0.8
        const totalNumBuckets = 16
        const bfSize = helper.computeBloomFilterSize(numElementsInserted, containmentScoreThresh, totalNumBuckets)
        expect(bfSize).toBeDefined()
    })

    it('returns a Bloom filter of optimal size with minimizers inserted', async () => {
        const seq = await testSeq.getSequence('1')
        const seqMinimizers = helper.extractMinimizers(seq)
        console.log('num minimizers in testSeq:', seqMinimizers.length)
        const containmentScoreThresh = 0.8
        const totalNumBuckets = 16
        const seqBF = helper.makeMinimizersBloomFilter(seqMinimizers, containmentScoreThresh, totalNumBuckets)
        expect(seqBF).toBeDefined()
    })
})
