/* Write bigsi and bucket map to file
 */
const config = require('../bigsi.config.json')
const fs = require('fs')

function computeNumBuckets(seqLength) {
    const numBuckets = Math.ceil(seqLength/config.bucketSize)
    return numBuckets
}

/** Makes the bucket map to retrieve sequence intervals from bucket number
 */
function makeBucketMap(seqName, seqSize, currentBucketNum, numBuckets){
    const bucketSize = config.bucketSize
    const bucketMap = {} 
    for (let bucketNum=0; bucketNum < numBuckets; bucketNum++){

        const intervalStart = Math.max(bucketNum*config.bucketSize - config.bucketOverhang, 0)
        let intervalEnd = Math.min(intervalStart + config.bucketSize + 2*config.bucketOverhang, seqSize)

        if (intervalStart === 0) { // handle first bucket 
            intervalEnd -= config.bucketOverhang
        }
        
        const adjustedBucketNum = bucketNum + currentBucketNum
        bucketMap[adjustedBucketNum] = {
            refName: seqName,
            bucketStart: intervalStart,
            bucketEnd: intervalEnd,
        }
    }

    return bucketMap
}

/** 
 * fasta - indexFasta of entire genome/sequence
 */
async function makeBucketToPositionMap(fasta){
    const seqSizes = await fasta.getSequenceSizes()

    let bucketToPositionMap = {}
    let currentBucketNum = 0
    for (seqName in seqSizes) {
        const seqSize = seqSizes[seqName]
        const numBuckets = computeNumBuckets(seqSize)
        const bucketMap = makeBucketMap(seqName, seqSize, currentBucketNum, numBuckets)
        bucketToPositionMap = { ...bucketToPositionMap, ...bucketMap }
        currentBucketNum += numBuckets
    }

    return bucketToPositionMap
 }

function writeBucketMapToJSON(bucketMap, output){
    // convert JSON object to string
    const json = JSON.stringify(bucketMap);

    // write JSON string to a file
    fs.writeFile(output, json, (err) => {
        if (err) {
            throw err;
        }
        console.log("Saved JSON of bucket map.");
    });
}

function writeQueryConfigToJSON(bigsiDims, output) {
    // convert JSON object to string
    const json = JSON.stringify(bigsiDims);

    // write JSON string to a file
    fs.writeFile(output, json, (err) => {
        if (err) {
            throw err;
        }
        console.log("Saved JSON of query configuration.");
    });
}

// Write BIGSI to binary file
// --------------------------

/** Converts bigsi into a flat array of bitstrings
 * @param { matrix } bigsi 
 */
function bigsiToBitstrings(bigsi) {
    const bigsiArr = bigsi()

    const bigsiBitstrings = []
    for (const row of bigsiArr){
        const rowBitstring = row.join('')
        bigsiBitstrings.push(rowBitstring)
    }

    return bigsiBitstrings
}

function bitstringsToInts(bitstrings) {
    const ints = []
    const intSize = config.intBits
    for (const bitstring of bitstrings) {
        const paddedBitstringSize = bitstring.length + (intSize - (bitstring.length % intSize))
        const paddedBitstring = bitstring.padEnd(paddedBitstringSize, '0')
        const regex = new RegExp(`.{1,${intSize}}`, "g");
        const chunks = paddedBitstring.match(regex)
        for (const chunk of chunks) {
            const integer = parseInt(chunk, 2)
            ints.push(integer)
        }
    }
    return ints
}

/** Converts bigsi into a flat array of ints
 */
function bigsiToInts(bigsi){
    console.log(`bigsi dimensions: ${bigsi.size()} [rows, cols]`)
    const bigsiArr = bigsi()

    const bigsiInts = []
    for (const row of bigsiArr){
        for (let i=0; i < row.length; i+=config.intBits){
            const seqRow = row.slice(i, i+config.intBits)
            const bitstring = seqRow.join('')
            const seqInt = parseInt(bitstring, 2)
            bigsiInts.push(seqInt)
        }
    }

    console.log('total number of ints:', bigsiInts.length)

    return bigsiInts
}

function flattenBigsi(bigsi) {
    const bigsiArr = bigsi()
    return bigsiArr.flat()
}

/**
 * @param { array } bigsiInts - bigsi matrix as a flat array of 16-bit ints
 *
 * @returns { TypedArray } binaryBigsi - bigsi as a TypedArray
 */
function makeBinaryBigsi(bigsiInts){
    const totalBufSize = Math.ceil(bigsiInts.length*2) // buffer size in bytes
    console.log('buffer size: ', totalBufSize)
    const buf = new ArrayBuffer(totalBufSize); // specifies the number of bytes in the buffer
    const bits = new Uint16Array(buf); // sets up a binary view on the buffer
    
    bits.set(bigsiInts); // adds the elements to the buffer based on the view
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

module.exports = {
    makeBucketMap: makeBucketMap,
    makeBucketToPositionMap: makeBucketToPositionMap,
    writeBucketMapToJSON: writeBucketMapToJSON,
    writeQueryConfigToJSON: writeQueryConfigToJSON,
    bitstringsToInts: bitstringsToInts,
    bigsiToBitstrings: bigsiToBitstrings,
    bigsiToInts: bigsiToInts,
    flattenBigsi: flattenBigsi,
    writeBinaryBigsi: writeBinaryBigsi,
    makeBinaryBigsi: makeBinaryBigsi,
    makeHexBigsi: makeHexBigsi,
}
