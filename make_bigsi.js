const commonFunc = require('./common_func.js')
const matrix = require('matrix-js')
const BitSet = require('bitset')
const fs = require('fs')

/**
 * @param { number } seqSize - size of the sequence
 * @param { string } seqName - name of sequence
 * @param { number } bucketNumStart, bucketNumEnd - first and last bucket 
 * numbers of the sequence
 * @param { number } overhang - bucket overhangs (equal to upper bound of query 
 * sequence length
 *
 * @returns { object } bucketMap - maps buckets => { refName: "1", bucketStart: 0, bucketEnd: 0 }
 */
function makeBucketMap(seqSize, seqName, bucketNumStart, bucketNumEnd, overhang=190000){
    const numBuckets = bucketNumEnd - bucketNumStart + 1
    const bucketSize = Math.round((seqSize + ((numBuckets - 1)*overhang))/numBuckets)

    const bucketMap = {} 
    for (let bucketNum=bucketNumStart; bucketNum <= bucketNumEnd; bucketNum++){
        bucketMap[bucketNum] = {}

        const intervalStart = Math.max((bucketSize-overhang)*(bucketNum%numBuckets), 0)
        let intervalEnd = intervalStart + bucketSize

        // Handle last bucket
        if (bucketNum == numBuckets - 1){
            intervalEnd = seqSize
        }

        bucketMap[bucketNum] = Object.assign(
            {
                refName: seqName,
                bucketStart: intervalStart,
                bucketEnd: intervalEnd,
            }, 
            bucketMap[bucketNum]
        )
        
    }

    return bucketMap
}

/** 
 * @param { indexedFasta } seq - indexFasta of entire genome/sequence
 * @param { number } seqSizeThreshold - min size of sequence (to filter small 
 * sequences)
 * @param { number } numBuckets - number of buckets per sequence
 *
 * @returns { object } bucketToPosition - maps buckets => { refName: "1", start: 0, end: 0 }
 */
async function makeBucketToPosition(seq, seqSizeThreshold, numBuckets=16){
    const seqSizes = await seq.getSequenceSizes()

    let bucketToPosition = {}
    let seqIdx = 0
    for (seqName in seqSizes) {
        const seqSize = seqSizes[seqName]
        if ( seqSize > seqSizeThreshold ) {
            const bucketStart = seqIdx*numBuckets
            const bucketEnd = bucketStart + numBuckets - 1
            const bucketMap = makeBucketMap(seqSize, seqName, bucketStart, bucketEnd)
            bucketToPosition = { ...bucketToPosition, ...bucketMap }
            seqIdx++
        }
    }

    return bucketToPosition
 }

/** 
 * @param { object } bucketToPosition - maps buckets => { refName: "1", start: 0, end: 1 }
 * @param { string } outputFilename - name of output file
 *
 * @returns json file with bucketToPosition
 */
function writeBucketMapToJSON(bucketToPosition, outputFilename){
    
    const bucketMapJSON = JSON.stringify(bucketToPosition, null, 4);

    fs.writeFile(outputFilename, bucketMapJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`Bucket-Position map written to: ${outputFilename}`);
        }
    });

}

async function splitSeqIntoBuckets(seq, seqName, numBuckets, overhang=300_000){
    /* Inputs:
     *  seq -- sequence object from IndexedFasta lib
     *  seqName -- name of sequence in fai file
     *  numBuckets -- number of overlapping buckets to split sequence into
     *
     * Outputs:
     *  buckets -- array of sequence strings
    */

    const seqSize = await seq.getSequenceSize(seqName);
    const bucketSize = Math.round((seqSize + ((numBuckets - 1)*overhang))/numBuckets)

    const buckets = []
    for (let i=0; i < numBuckets; i++){
        const bucketStart = Math.max((bucketSize-overhang)*i, 0)
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

function makeBucketsBloomFilters(bucketSequences){
    /* Inputs:
     *  bucketSequences -- array of bucket sequences (strings)
     *
     * Outputs:
     *  bucketsBloomFilters -- array of Bloom filters (one for each bucket)
     *
     */

    const bucketsBloomFilters = []
    for (let idx=0; idx < bucketSequences.length; idx++){
        const bucketSequence = bucketSequences[idx]
        const bucketMinimizers = commonFunc.extractMinimizers(bucketSequence)

        const bucketBloomFilter = commonFunc.makeMinimizersBloomFilter(bucketMinimizers)
        console.log(`${bucketBloomFilter.length} minimizers inserted into bloom filter for bucket ${idx}`)
        bucketsBloomFilters.push(bucketBloomFilter)
    }

    console.log('total number of bloom filters: ', bucketsBloomFilters.length)
    return bucketsBloomFilters
}

async function buildBigsi(bucketSequences){
    /* Inputs:
     *  bucketSequences -- array of bucket sequences (strings) 
     *
     * Outputs:
     *  bigsi matrix of bucket Bloom filters
     *
     */

    const bucketsBloomFilters = makeBucketsBloomFilters(bucketSequences)
    const bigsiMatrix = matrix(bucketsBloomFilters.map(bloomFilterObj => bloomFilterObj._filter))
    const bigsi = matrix(bigsiMatrix.trans())

    bucketSequences = []

    return bigsi

}

/**
 * @param { String } genomesDir - directory where genomes are stored (in fasta 
 * format with fai indices)
 *
 * @returns { array of strings } genomeSeqs -- array of genome sequences 
 * (each sequence is a bucket)
 */
async function getMultiGenomeSeqs(genomesDir){
    const dir = await fs.promises.opendir(genomesDir)

    const genomeSeqs = []
    for await (const file of dir) {
        const isFasta = file.name.endsWith('.fasta') || file.name.endsWith('.fa') || file.name.endsWith('.fsa')

        if (isFasta){
            console.log('Opening:', file.name)
            const fai = file.name + '.fai'
            const genome = await commonFunc.loadFasta(`${dir.path}/${file.name}`, `${dir.path}/${fai}`)
            const seqSizeThreshold = 0
            const seqNames = await commonFunc.getFilteredGenomeSeqs(genome, seqSizeThreshold)

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
 * @param { IndexedFasta } genome - indexedfasta object of genome
 *
 * @returns { array of matrices } genomeBigsis - array of bigsi matrices for 
 *  each seq in genome
 */
async function makeGenomeBigsis(genome, numBuckets, seqSizeThreshold=10**7){
    const seqNames = await commonFunc.getFilteredGenomeSeqs(genome, seqSizeThreshold)
    console.log('seqNames: ', seqNames)

    const genomeBigsis = []
    for (let i=0; i < seqNames.length; i++){
        let seqBuckets = await splitSeqIntoBuckets(genome, seqNames[i], numBuckets) 
        console.log(`Split ${seqNames[i]} into buckets, building bigsi...`)
        const seqBigsi = await buildBigsi(seqBuckets)
        console.log(`Bigsi of ${seqNames[i]} built, pushing into array...`)
        genomeBigsis.push(seqBigsi)
        const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
        console.log(`Process uses ${memoryUsed} MB`)
        seqBuckets = []
    }

    return genomeBigsis
}

/**
 * @param { array of functions } bigsis - array of bigis matrices
 *
 * @returns { function } mergedBigsi - merged bigsi
 */

async function mergeBigsis(bigsis){
    let mergedBigsi = bigsis[0]
    for (let i=1; i < bigsis.length; i++){
        mergedBigsi = matrix(mergedBigsi.merge.right(bigsis[i]()))
    }

    return mergedBigsi
}


/** Converts bigsi into a flat array of ints
 */
function bigsiToInts(bigsi, intSize){
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
 * @returns { TypedArray } binaryDumpBigsi - bigsi as a TypedArray
 */
function makeBinaryDumpBigsi(u16IntRows){
    const totalBufSize = Math.ceil(u16IntRows.length*2) // buffer size in bytes
    console.log('buffer size: ', totalBufSize)
    const buf = new ArrayBuffer(totalBufSize); // specifies the number of bytes in the buffer
    const bits = new Uint16Array(buf); // sets up a binary view on the buffer
    
    bits.set(u16IntRows); // adds the elements to the buffer based on the view
    console.log('num elements in u16IntArray: ', bits.length)
    
    return bits
}



/**
 * @param { TypedArray } binaryDumpBigsi - bigsi as a TypedArray 
 * @param { string } filename - name of the file we are writing to
 *
 * @returns - writes TypedArray to binary file
 */
function writeBinaryDumpBigsi(binaryDumpBigsi, filename){

    fs.writeFileSync(filename, binaryDumpBigsi, 'binary', function (err) {
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI binary dump written to: ${filename}`);
            }
    });
}

/**
 * @param { Array of Arrays } bigsi
 *
 * @returns { Array of Hexdecimal strings } hexBigsi
 */
function makeHexBigsi(bigsi){
    let hexBigsi = []
    const numRows = bigsi.size()[0]

    for (let i = 0; i < numRows; i++){
        const ithRow = bigsi(i)
        const binaryString = ithRow.join('')
        const bitset = BitSet.fromBinaryString(binaryString)
        const hexString = bitset.toString(16)
        hexBigsi.push(hexString)
    }

    return hexBigsi
}


/** 
 * @param { Array } bigsi -- bigsi matrix object as array of strings
 * @param { string } outputFilename -- output file name string
 *
 * @returns writes bigsi to a json file
 */
function writeBigsiToJSON(bigsi, outputFilename){

    const bigsiJSON = JSON.stringify(bigsi, null, 4);

    fs.writeFile(outputFilename, bigsiJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI written to: ${outputFilename}`);
        }
    });

}


module.exports = {
    makeBucketMap: makeBucketMap,
    makeBucketToPosition: makeBucketToPosition,
    writeBucketMapToJSON: writeBucketMapToJSON,
    splitSeqIntoBuckets: splitSeqIntoBuckets,
    makeBucketsBloomFilters: makeBucketsBloomFilters,
    getMultiGenomeSeqs: getMultiGenomeSeqs,
    buildBigsi: buildBigsi,
    makeGenomeBigsis: makeGenomeBigsis,
    mergeBigsis: mergeBigsis,
    bigsiToInts: bigsiToInts,
    makeBinaryDumpBigsi: makeBinaryDumpBigsi,
    writeBinaryDumpBigsi: writeBinaryDumpBigsi,
    makeHexBigsi: makeHexBigsi,
    writeBigsiToJSON: writeBigsiToJSON,
    makeBucketToGenome: makeBucketToGenome
}
