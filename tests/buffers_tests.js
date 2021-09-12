const hexBigsi = require('./data/hg38_chr1and2_hexbigsi.json')
const BitSet = require('bitset')

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

function main(){
    intToBitstringTest()
}

main()
