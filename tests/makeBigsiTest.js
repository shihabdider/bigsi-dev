const makeBigsi = require('../make_bigsi.js')
const { IndexedFasta } = require('@gmod/indexedfasta')

// Init a seq object (hg38 chr1 ~3MB subseq) for use with tests
const fastaPath = './data/testRefSeq.fa'
const faiPath = './data/testQuerySeq.fa.fai'

const SEQ = new IndexedFasta({
    path: fastaPath,
    faiPath: faiPath,
    chunkSizeLimit: 50000000
});

async function makeBigsiTest(){
    const seqSize = await SEQ.getSequenceSize('1')
    console.log(seqSize)
    const bucketSequences = await makeBigsi.splitSeqIntoBuckets(SEQ, seqName='1', overhang=3000, numBuckets=5)
    const bigsiMatrix = await makeBigsi.buildBigsi(bucketSequences)
    console.log(bigsiMatrix(1))
}

async function main(){
    await makeBigsiTest()
}

main()
