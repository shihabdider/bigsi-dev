const Uint1Array = require('uint1array').default;
const BitSet = require('bitset')
const matrix = require('matrix-js')


const NUM_ROWS = 160_000

/**
 * @param {number} numCols - number of buckets/columns in bigsi
 * @return {array of arrays} - bigsi as an array of arrays
 */
function makeRandomBigsi(numCols){
    bigsi = []
    for (let i = 0; i < NUM_ROWS; i++){
        const bs = BitSet.Random(numCols);
        //bigsi.push(bs.toString().padStart(numCols, '0'))
        bigsi.push(bs.toString().padStart(numCols, '0').split('').map(Number))
    }

    return bigsi
}

/** Converts bigsi into a flat array of u16ints
 */
function bigsiToInts(bigsi, intSize){
    const flatRows = bigsi().flat()
    console.log('length and num ints:', flatRows.length, flatRows.length/intSize)

    const intRows = []
    for (let start=0; start < flatRows.length/intSize; start+=intSize){
        const end = start + intSize
        const row = flatRows.slice(start, end)
        const bitstring = row.join('')
        const rowInt = parseInt(bitstring, 2)
        intRows.push(rowInt)
    }

    return intRows
}

function bigsiTypedArray(u16IntRows, numCols){
    const totalBufSize = Math.ceil(numCols*NUM_ROWS/8) 
    const buf = new ArrayBuffer(totalBufSize); // specifies the number of bytes in the buffer
    const bits = new Uint16Array(buf); // sets up a binary view on the buffer
    
    bits.set(u16IntRows); // adds the elements to the buffer based on the view
    
    return bits
}

function queryBinaryBigsi(binaryBigsi, rowNum){
    const offsetStart = rowNum*2
    const offsetEnd = offsetStart + 2
    console.log(`From binary bigsi: ${binaryBigsi.slice(offsetStart, offsetEnd).join(',')}`, binaryBigsi.length);

    const rowBitString = Array.from(binaryBigsi.slice(offsetStart, offsetEnd)).map( (num) => {return num.toString(2)}).join('')
    console.log('binary dump row as bitstring: ', rowBitString)
    const bs = new BitSet(rowBitString)
    console.log('bitset: ', bs.toString())
}

function main(){
    const numCols = 32
    const intSize = 16
    const randBigsi = matrix(makeRandomBigsi(numCols))
    const u16IntRows = bigsiToInts(randBigsi, intSize)
    //const mergedRows = randBigsi().flat()
    
    const bits = bigsiTypedArray(u16IntRows, numCols)
    const randRow = Math.floor(Math.random()*16000)
    //const randRow = 1
    queryBinaryBigsi(bits, randRow)
    const rawBigsiRow = [
        parseInt(randBigsi(randRow).slice(0,16).join(''), 2), 
        parseInt(randBigsi(randRow).slice(16,32).join(''), 2)]
    console.log('from raw Bigsi: ', rawBigsiRow)
    console.log('from raw Bigsi (binary): ', randBigsi(randRow))
}

main()
