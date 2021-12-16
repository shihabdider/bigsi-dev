/* Write bigsi and bucket map to file
 */
const config = require('../bigsi.config.json')
const helper = require('./helper.js')
const fs = require('fs')

function makeBucketMap(seqName, seqSize, seqIdx){
    const bucketStart = seqIdx*config.numBuckets
    const bucketEnd = bucketStart + config.numBuckets - 1
    const bucketSize = Math.round(
        (seqSize + ((config.numBuckets - 1)*config.bucketOverhang))/config.numBuckets
    )

    const bucketMap = {} 
    for (let bucketNum=bucketStart; bucketNum <= bucketEnd; bucketNum++){

        const intervalStart = Math.max(
            (bucketSize-config.bucketOverhang)*(bucketNum%config.numBuckets), 0
        )
        let intervalEnd = intervalStart + bucketSize

        // Handle final bucket
        if (bucketNum == config.numBuckets - 1){
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
 * fasta - indexFasta of entire genome/sequence
 */
async function makeBucketToPositionMap(fasta){
    const seqSizes = await fasta.getSequenceSizes()

    let bucketToPositionMap = {}
    let seqIdx = 0
    for (seqName in seqSizes) {
        const seqSize = seqSizes[seqName]
        const bucketMap = makeBucketMap(seqName, seqSize, seqIdx)
        bucketToPositionMap = { ...bucketToPositionMap, ...bucketMap }
        seqIdx++
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

// Write BIGSI to binary file
// --------------------------

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
    bigsiToInts: bigsiToInts,
    writeBinaryBigsi: writeBinaryBigsi,
    makeBinaryBigsi: makeBinaryBigsi,
    makeHexBigsi: makeHexBigsi,
}
