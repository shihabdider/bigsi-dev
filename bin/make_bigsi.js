/* Makes a bigsi from a reference sequence
 *
 * Input:
 *  sequence: (string)
 *  bloom filter size: (number)
 *  number of buckets: (number)
 *  size of buckets: (number)
 *
 * Output:
 *  bigsi: (matrix)
 *
 */
const helper = require('./helper.js')
const matrix = require('matrix-js')
const BitSet = require('bitset')
const fs = require('fs')

// Globals
const BUCKET_OVERHANG = 150_000 //half of max query sequence length

function makeBucketMap(seqName, seqSize, seqIdx, numBucketsPerSeq){
    const bucketStart = seqIdx*numBucketsPerSeq
    const bucketEnd = bucketStart + numBucketsPerSeq - 1
    const bucketSize = Math.round((seqSize + ((numBucketsPerSeq - 1)*BUCKET_OVERHANG))/numBucketsPerSeq)

    const bucketMap = {} 
    for (let bucketNum=bucketStart; bucketNum <= bucketEnd; bucketNum++){

        const intervalStart = Math.max((bucketSize-BUCKET_OVERHANG)*(bucketNum%numBucketsPerSeq), 0)
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

async function splitSeqIntoBuckets(seq, seqName, numBuckets){

    const seqSize = await seq.getSequenceSize(seqName);
    const bucketSize = Math.round((seqSize + ((numBuckets - 1)*BUCKET_OVERHANG))/numBuckets)

    const buckets = []
    for (let i=0; i < numBuckets; i++){
        const bucketStart = Math.max((bucketSize-BUCKET_OVERHANG)*i, 0)
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

function makeBucketBloomFilter(sequence, bloomFilterSize){
    const bucketMinimizers = helper.extractMinimizers(sequence)
    const bucketBloomFilter = helper.makeMinimizersBloomFilter(
            bucketMinimizers, 
            bloomFilterSize
        )

    return bucketBloomFilter
}

async function buildBigsi(sequence, bloomFilterSize, bucketCoords){

    const seqBloomFilters = []
    for (const coord in bucketCoords){
        const ithBucketSequence = sequence.slice(coord.bucketStart, coord.bucketEnd);
        console.log(coord.bucketStart, coord.bucketEnd, ithBucketSequence.length)
        const bucketBloomFilter = makeBucketBloomFilter(ithBucketSequence, bloomFilterSize)
    }

    let bigsi = matrix(seqBloomFilters.map(bloomFilterObj => bloomFilterObj._filter))
    bigsi = matrix(bigsi.trans()) // transpose to make bloom filters into columns of matrix

    return bigsi

}

function computeNumMinimizers(seqLength, windowSize=100){
    const numMinimizers = Math.ceil(seqLength/windowSize * 2)
    return numMinimizers
}

function computeBloomFilterSize(maxNumElementsInserted, containmentScoreThresh, totalNumBuckets){
    // initialize set parameters
    const minQueryMinimizers = 100  // 5Kbp min query = 100 minimizers
    const falseHitThresh = 1e-2
    // iterate over a array size range...
    for ( let bloomFilterSize = 0; bloomFilterSize <= 1e6; bloomFilterSize += 1e3 ){
        const falsePosRate = computeBloomFilterFalsePosRate(maxNumElementsInserted, bloomFilterSize)
        const falseHitProb = computeFalseHitProb(
            falsePosRate, 
            minQueryMinimizers, 
            containmentScoreThresh
        )

        // accounting for all buckets in bigsi
        const falseHitProbUpper = falseHitProb*totalNumBuckets
        // break if false hit rate less than threshold and return
        if ( falseHitProbUpper <= falseHitThresh ) {
            console.log(`optimal bloom filter size: ${bloomFilterSize}`)
            return bloomFilterSize
        }
    }
}

/* estimate bloom filter size using minimizer count computed from bucket size 
 * of longest sequence + bucket overhang
*/ 
function estimateBloomFilterSize(seqSizes, numBuckets){
    const seqSizesArr = Object.values(seqSizes)
    const maxSeqLength = Math.max(...seqSizesArr)/numBuckets + BUCKET_OVERHANG
    console.log('maximum sequence length: ', maxSeqLength)
    const maxNumElementsInserted = computeNumMinimizers(maxSeqLength)
    const containmentScoreThresh = 0.80
    const totalNumBuckets = seqSizes.length*numBuckets

    const bloomFilterSize = computeBloomFilterSize(
        maxNumElementsInserted, 
        containmentScoreThresh, 
        totalNumBuckets
    )

    return bloomFilterSize
}

function computeBucketCoords(seqLength, bucketSize, numBuckets) {
    const bucketCoords = []
    for (let i=0; i < numBuckets; i++){
        const bucketStart = Math.max((bucketSize-BUCKET_OVERHANG)*i, 0)
        let bucketEnd = bucketStart + bucketSize

        // Handle last bucket
        if (i == numBuckets - 1){
            bucketEnd = seqSize
        }

        const coord = { bucketStart, bucketEnd }
        bucketCoords.push(coord)
    }

    return bucketCoords
}

/**
 * @param { IndexedFasta } fasta - indexedFasta object
 *
 * @returns { matrix[] } fastaBigsis - bigsi matrix for each seq in fasta
 */
async function makeFastaBigsis(fasta, numBuckets){
    const seqNames = await fasta.getSequenceList()
    console.log('seqNames: ', seqNames)

    const seqSizes = await fasta.getSequenceSizes()
    const bloomFilterSize = estimateBloomFilterSize(seqSizes, numBuckets)

    const fastaBigsis = []
    for (const seqName of seqNames){
        const sequence = await fasta.getSequence(seqName)
        const bucketSize = Math.round(
            (seqSizes.seqName + ((numBuckets - 1)*BUCKET_OVERHANG))/numBuckets
        )
        const bucketCoords = computeBucketCoords(sequence.length, bucketSize, numBuckets)
        const seqBigsi = await buildBigsi(sequence, bloomFilterSize, bucketCoords)
        console.log(`Bigsi of ${seq} built...`)
        fastaBigsis.push(seqBigsi)
        const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
        console.log(`Process used ${memoryUsed} MB`)
    }

    return fastaBigsis
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

async function main(fasta, numBuckets) {
    const minSeqLength = 30e6
    const seqSizes = Object.values(await fasta.getSequenceSizes())
    const areFastaSeqsValidSize = Math.min(...seqSizes) > minSeqLength 

    if (areFastaSeqsValidSize && numBuckets > 0) {
        const bigsis = await makeFastaBigsis(fasta, numBuckets)
        console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
        const bigsi = await mergeBigsis(bigsis)
        console.log(`Bigsis merged!`)

        const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        return bigsi
        //if (isHexBigsi) {
        //    const hexBigsi = makeHexBigsi(bigsi)
        //    return hexBigsi
        //} else {
        //    const u16IntRows = bigsiToInts(bigsi)
        //    const binaryBigsi = makeBinaryBigsi(u16IntRows)

        //    return binaryBigsi
        //}
    } else { 
        if (!areFastaSeqsValidSize) { 
            console.log('All sequences must be at least 30Mbp in length.') 
        }
        if (!(numBuckets > 0)) { 
            console.log('Number of buckets must be greater than 0.') 
        }
    }
}

module.exports = {
    makeBucketMap: makeBucketMap,
    makeBucketToPositionMap: makeBucketToPositionMap,
    writeBucketMapToJSON: writeBucketMapToJSON,
    splitSeqIntoBuckets: splitSeqIntoBuckets,
    makeBucketBloomFilter: makeBucketBloomFilter,
    getMultiGenomeSeqs: getMultiGenomeSeqs,
    buildBigsi: buildBigsi,
    makeFastaBigsis: makeFastaBigsis,
    mergeBigsis: mergeBigsis,
    bigsiToInts: bigsiToInts,
    makeBinaryBigsi: makeBinaryBigsi,
    writeBinaryBigsi: writeBinaryBigsi,
    makeHexBigsi: makeHexBigsi,
    writeBigsiToJSON: writeBigsiToJSON,
    makeBucketToGenome: makeBucketToGenome,
    main: main
}
