const hexBigsi = require('./test_data/hg38_chr1_hex.json')
const helper = require('../bin/helper.js')
const BitSet = require('bitset')
const matrix = require('matrix-js')
const fs = require('fs')



//write bitmatrix to buffer 
function matrixToBuffer(bitmatrix){
    console.log(bitmatrix())
    const intSize = 16
    const flatRows = bitmatrix().flat()

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

    const totalBufSize = Math.ceil(intRows.length*2) // buffer size in bytes
    console.log('buffer size: ', totalBufSize)
    const buf = new ArrayBuffer(totalBufSize); // specifies the number of bytes in the buffer
    const bits = new Uint16Array(buf); // sets up a binary view on the buffer
    
    bits.set(intRows); // adds the elements to the buffer based on the view
    console.log('num elements in u16IntArray: ', bits.length)
    
    return bits
}

//load buffer file into bitmatrix
function bufferToMatrix(buffer){

    const matrixRows = []
    const numSeqs = 1
    const numRows = 2

    for (let i=0; i < numRows; i++){

        const offsetStart = i*numSeqs
        const offsetEnd = offsetStart + numSeqs
        console.log(offsetStart, offsetEnd)

        const rowInts = Array.from(buffer.subarray(offsetStart, offsetEnd))
        console.log(rowInts, rowInts[0].toString(2))
        const rowBitStrings = rowInts.map((num) => helper.zeroPadBitstring(num.toString(2), 16))
        const rowBitString = rowBitStrings.join('')

        // Front padding ensures all columns are accounted for
        const row = rowBitString.split('').map(Number)
        matrixRows.push(row)
    }

    const bitmatrix = matrix(matrixRows)

    return bitmatrix
}

function zeroPad(num, places){
    const paddedString = String(num).padStart(places, '0')
    return paddedString
}

function bitStringToUint8Array(bitString) {

    // split the string into octets
    const pairs = bitString.match(/[\d]{1,8}/gi)

    // convert the octets to integers
    const integers = pairs.map(function(s) {
        return parseInt(s, 2);
    });

    const array = new Uint8Array(integers);

    console.log(pairs, integers, array.buffer, array.length)
    return array;
}

function intToBitstringTest(){
    //const rowBitstring = '01000000000000000000001000000' //length = 20

    const toBinString = (bytes) => 
        bytes.reduce((str, byte) => str + byte.toString(2).padStart(8, '0'), '');
    
    for (const row of hexBigsi.slice(0,1)){
        const bs = new BitSet(`0x${row}`)
        const rowBitstring =  bs.toString()
        const Uint8Arr = bitStringToUint8Array(rowBitstring)
        console.log(rowBitstring, toBinString(Uint8Arr))
    }
    //console.log(toBinString(buffer))
    //console.log(rowBitstring, rowInt.toString(2))
}

function bufferTest(){
    const bufferSize = 3*300_000
    let buffer = new ArrayBuffer(bufferSize)
    let int8View = new Int8Array(buffer);

    console.log(int8View.length)
}

function writeBufferToFile(path, buffer){

}

function main(){
    const testBigsiArray = [
        [1, 1, 1, 0, 1, 1, 0, 1, 
            1, 1, 0, 1, 1, 0, 1, 1],
        [0, 0, 1, 0, 1, 1, 0, 1, 
            1, 1, 0, 0, 1, 0, 1, 0]
    ]

    const bigsiMatrix = matrix(testBigsiArray)
    const buffer = matrixToBuffer(bigsiMatrix)

    const bufferPath ='tests/test_data/bufferTestFile.bin' 
    fs.writeFileSync(bufferPath, buffer, 'binary')
    const loadedBuffer = fs.readFileSync(bufferPath)
    let arr = new Uint16Array(loadedBuffer.buffer, loadedBuffer.byteOffset, loadedBuffer.length / 2);
    console.log('arr', arr)

    const derivedMatrix = bufferToMatrix(arr)

    console.log('bigsi:', bigsiMatrix(), '\n\n', 'derived:', derivedMatrix(), '\n\n', 'original:', testBigsiArray)

}

main()
