const helper = require('./helper.js')
const matrix = require('matrix-js')
const BitSet = require('bitset')
const fs = require('fs')

// bucketOverhang - overhangs (equal to half the upper bound of query sequence length or 150K) 
function makeBucketMap(seqName, seqSize, seqIdx, numBucketsPerSeq, bucketOverhang=150_000){
    const bucketStart = seqIdx*numBucketsPerSeq
    const bucketEnd = bucketStart + numBucketsPerSeq - 1
    const bucketSize = Math.round((seqSize + ((numBucketsPerSeq - 1)*bucketOverhang))/numBucketsPerSeq)

    const bucketMap = {} 
    for (let bucketNum=bucketStart; bucketNum <= bucketEnd; bucketNum++){

        const intervalStart = Math.max((bucketSize-bucketOverhang)*(bucketNum%numBucketsPerSeq), 0)
        let intervalEnd = intervalStart + bucketSize

        // Handle final bucket
        if (bucketNum == numBucketsPerSeq - 1){
            intervalEnd = seqSize
        }

        bucketMap[bucketNum] = {
                refName: seqName,
                bucketStart: intervalStart,
                bucketEnd: intervalEnd,
            }
    }

    return bucketMap
}

/** 
 * seq - indexFasta of entire genome/sequence
 * seqSizeThreshold - min size of sequence (to filter small sequences)
 */
async function makeBucketToPositionMap(seq, seqSizeThreshold, numBucketsPerSeq=16){
    const seqSizes = await seq.getSequenceSizes()

    let bucketToPositionMap = {}
    let seqIdx = 0
    for (seqName in seqSizes) {
        const seqSize = seqSizes[seqName]
        if ( seqSize >= seqSizeThreshold ) {
            const bucketMap = makeBucketMap(seqName, seqSize, seqIdx, numBucketsPerSeq)
            bucketToPositionMap = { ...bucketToPositionMap, ...bucketMap }
            seqIdx++
        }
    }

    return bucketToPositionMap
 }

function writeBucketMapToJSON(bucketToPositionMap, outputFilename){
    
    const bucketMapJSON = JSON.stringify(bucketToPositionMap, null, 4);

    fs.writeFile(outputFilename, bucketMapJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`Bucket-Position map written to: ${outputFilename}`);
        }
    });

}
/**
 * @param { String } genomesDir - directory where genomes are stored (in fasta 
 * format with fai indices)
 *
 * @returns { Object } bucketToGenome -- maps bucket/column of bigsi to 
 * corresponding genome filename
 */
async function makeBucketToGenome(genomesDir){
    const dir = await fs.promises.opendir(genomesDir)

    const bucketToGenome = {}
    const currentCol = 0
    for await (const file of dir) {
        const isFasta = file.name.endsWith('.fasta') || file.name.endsWith('.fa')
        if (isFasta){
            const genomeName = file.name.replace(/\.[^/.]+$/, "")
            bucketToGenome[currentCol] = genomeName
            currentCol++
        }
    }

    return bucketToGenome
}

/**
 * @param { String } genomesDir - directory where genomes are stored (in fasta 
 * format with fai indices)
 *
 * @returns { array of strings } genomeSeqs -- array of genome sequences 
 * (each sequence corresponds to a bucket)
 */
async function getMultiGenomeSeqs(genomesDir){
    const dir = await fs.promises.opendir(genomesDir)

    const genomeSeqs = []
    for await (const file of dir) {
        const isFasta = file.name.endsWith('.fasta') || file.name.endsWith('.fa') || file.name.endsWith('.fsa')

        if (isFasta){
            console.log('Opening:', file.name)
            const fai = file.name + '.fai'
            const genome = await helper.loadFasta(`${dir.path}/${file.name}`, `${dir.path}/${fai}`)
            const seqSizeThreshold = 0
            const seqNames = await helper.getFilteredGenomeSeqs(genome, seqSizeThreshold)

            const seqStrings = [] 
            for (const seqName of seqNames){
                seq = await genome.getSequence(seqName);   
                seqStrings.push(seq)
            }
            const genomeStr = seqStrings.join('')
            console.log(`length of genome ${file.name}:`, genomeStr.length)
            genomeSeqs.push(genomeStr)
        }
    }
    
    return genomeSeqs
}

async function splitSeqIntoBuckets(seq, seqName, numBuckets, bucketOverhang=150_000){

    const seqSize = await seq.getSequenceSize(seqName);
    const bucketSize = Math.round((seqSize + ((numBuckets - 1)*bucketOverhang))/numBuckets)

    const buckets = []
    for (let i=0; i < numBuckets; i++){
        const bucketStart = Math.max((bucketSize-bucketOverhang)*i, 0)
        let bucketEnd = bucketStart + bucketSize

        // Handle last bucket
        if (i == numBuckets - 1){
            bucketEnd = seqSize
        }

        const nthBucketSequence = await seq.getSequence(seqName, bucketStart, bucketEnd);
        buckets.push(nthBucketSequence)
        console.log(bucketStart, bucketEnd, nthBucketSequence.length)
    
    }
    return buckets
}

function makeBucketsBloomFilters(maxNumElementsInserted, bucketSequences, totalNumBuckets){

    const bucketsMinimizers = []

    console.log(`Max num elements inserted: ${maxNumElementsInserted}`)

    const containmentScoreThresh = 0.80
    const bloomFilterSize = helper.computeBloomFilterSize(
        maxNumElementsInserted, 
        containmentScoreThresh, 
        totalNumBuckets
    )

    const bucketsBloomFilters = []
    for (bucketMinimizers of bucketsMinimizers){
        const bucketBloomFilter = helper.makeMinimizersBloomFilter(
                bucketMinimizers, 
                bloomFilterSize
            )
        bucketsBloomFilters.push(bucketBloomFilter)
    }

    console.log('total number of bloom filters: ', bucketsBloomFilters.length)
    return bucketsBloomFilters
}

async function buildBigsi(bucketSequences, maxNumElementsInserted){

    const totalNumBuckets = bucketSequences.length
    const bucketsBloomFilters = makeBucketsBloomFilters(bucketSequences, totalNumBuckets, maxNumElementsInserted)
    const bigsiMatrix = matrix(bucketsBloomFilters.map(bloomFilterObj => bloomFilterObj._filter))
    const bigsi = matrix(bigsiMatrix.trans())

    return bigsi

}

/**
 * @param { IndexedFasta } genome - indexedFasta object of genome
 *
 * @returns { array of matrices } genomeBigsis - array of bigsi matrices for 
 *  each seq in genome
 */
async function makeGenomeBigsis(genome, numBuckets, seqSizeThreshold=3*10**7){
    const seqNames = await helper.getFilteredGenomeSeqs(genome, seqSizeThreshold)
    console.log('seqNames: ', seqNames)

    // use estimate of minimizer count (for computing bloom filter size) based 
    // on longest sequence
    const seqSizes = await Object.values(genome.getSequenceSizes())
    const maxSeqLength = Math.max(seqSizes)
    const maxNumElementsInserted = helper.computeNumMinimizers(maxSeqLength)

    const genomeBigsis = []
    for (let i=0; i < seqNames.length; i++){
        let seqBuckets = await splitSeqIntoBuckets(genome, seqNames[i], numBuckets) 
        console.log(`Split ${seqNames[i]} into buckets, building bigsi...`)
        const seqBigsi = await buildBigsi(seqBuckets, maxNumElementsInserted)
        console.log(`Bigsi of ${seqNames[i]} built, pushing into array...`)
        genomeBigsis.push(seqBigsi)
        const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
        console.log(`Process uses ${memoryUsed} MB`)
        seqBuckets = null
    }

    return genomeBigsis
}

function mergeBigsis(bigsis){
    let mergedBigsi = bigsis[0]
    if (bigsis.length > 1){
        for (let i=1; i < bigsis.length; i++){
            mergedBigsi = matrix(mergedBigsi.merge.right(bigsis[i]()))
        }
    }

    return mergedBigsi
}


/** Converts bigsi into a flat array of ints
 */
function bigsiToInts(bigsi, intSize=16){
    console.log(`bigsi dimensions: ${bigsi.size()} [rows, cols]`)
    const flatRows = bigsi().flat()

    const intRows = []
    for (let i=0; i < flatRows.length/intSize; i++){
        const start = i*intSize
        const end = start + intSize
        const row = flatRows.slice(start, end)
        const bitstring = row.join('')
        const rowInt = parseInt(bitstring, 2)
        intRows.push(rowInt)
    }

    console.log('flat bigsi length and num ints:', flatRows.length, intRows.length)

    return intRows
}


/**
 * @param { array } u16IntRows - bigsi matrix as a flat array of 16-bit ints
 *
 * @returns { TypedArray } binaryBigsi - bigsi as a TypedArray
 */
function makeBinaryBigsi(u16IntRows){
    const totalBufSize = Math.ceil(u16IntRows.length*2) // buffer size in bytes
    console.log('buffer size: ', totalBufSize)
    const buf = new ArrayBuffer(totalBufSize); // specifies the number of bytes in the buffer
    const bits = new Uint16Array(buf); // sets up a binary view on the buffer
    
    bits.set(u16IntRows); // adds the elements to the buffer based on the view
    console.log('num elements in u16IntArray: ', bits.length)
    
    return bits
}



/** 
 * @param { TypedArray } binaryBigsi - bigsi matrix in TypedArray format
 * @param { string } outputPath
 */
function writeBinaryBigsi(binaryBigsi, outputPath){

    fs.writeFileSync(outputPath, binaryBigsi, 'binary', function (err) {
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI binary dump written to: ${outputPath}`);
            }
    });
}

function makeHexBigsi(bigsi){
    const numRows = bigsi.size()[0]

    let hexBigsi = []
    for (let i = 0; i < numRows; i++){
        const ithRow = bigsi(i)
        const binaryString = ithRow.join('')
        const bitset = BitSet.fromBinaryString(binaryString)
        const hexString = bitset.toString(16)
        hexBigsi.push(hexString)
    }

    return hexBigsi
}


function writeBigsiToJSON(bigsi, outputPath){

    const bigsiJSON = JSON.stringify(bigsi, null, 4);

    fs.writeFile(outputPath, bigsiJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI written to: ${outputPath}`);
        }
    });

}

async function main(seq, numBuckets, seqSizeThreshold, isHexBigsi=false){
    const bigsis = await makeGenomeBigsis(seq, numBuckets, seqSizeThreshold)
    console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
    const bigsi = await mergeBigsis(bigsis)
    console.log(`Bigsis merged!`)

    const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
    console.log(`Process uses ${memoryUsed}`)

    if (isHexBigsi) {
        const hexBigsi = makeHexBigsi(bigsi)
        return hexBigsi
    } else {
        const u16IntRows = bigsiToInts(bigsi)
        const binaryBigsi = makeBinaryBigsi(u16IntRows)

        return binaryBigsi
    }
}

module.exports = {
    makeBucketMap: makeBucketMap,
    makeBucketToPositionMap: makeBucketToPositionMap,
    writeBucketMapToJSON: writeBucketMapToJSON,
    splitSeqIntoBuckets: splitSeqIntoBuckets,
    makeBucketsBloomFilters: makeBucketsBloomFilters,
    getMultiGenomeSeqs: getMultiGenomeSeqs,
    buildBigsi: buildBigsi,
    makeGenomeBigsis: makeGenomeBigsis,
    mergeBigsis: mergeBigsis,
    bigsiToInts: bigsiToInts,
    makeBinaryBigsi: makeBinaryBigsi,
    writeBinaryBigsi: writeBinaryBigsi,
    makeHexBigsi: makeHexBigsi,
    writeBigsiToJSON: writeBigsiToJSON,
    makeBucketToGenome: makeBucketToGenome,
    main: main
}
