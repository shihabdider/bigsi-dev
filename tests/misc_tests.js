const BitSet = require('bitset')
const fs = require('fs')
const Uint1Array = require('uint1array').default;
const commonFunc = require('../common_func.js')
const matrix = require('matrix-js')


function bitStringToArrayBuffer(bitString) {

    // split the string into octets
    var pairs = bitString.match(/[\d]{1,8}/gi)

    // convert the octets to integers
    var integers = pairs.map(function(s) {
        return parseInt(s, 2);
    });

    var array = new Uint8Array(integers);

    return array.buffer;
}

/**
 * Convert a hex string to an ArrayBuffer.
 * 
 * @param {string} hexString - hex representation of bytes
 * @return {ArrayBuffer} - The bytes in an ArrayBuffer.
 */
function hexStringToArrayBuffer(hexString) {
    // remove the leading 0x
    hexString = hexString.replace(/^0x/, '');
    
    // split the string into pairs of octets
    var pairs = hexString.match(/[\dA-F]{2}/gi);
    
    // convert the octets to integers
    var integers = pairs.map(function(s) {
        return parseInt(s, 16);
    });
    
    var array = new Uint8Array(integers);
    
    return array.buffer;
}

const zeroPad = (num, places) => String(num).padStart(places, '0')

function makeBigsiRowArrBuf(){

    bigsi = []
    for (let i = 0; i < 300000; i++){
        const bs = BitSet.Random(230);
        bigsi.push(bs.toString())
    }

    for (let i = 0; i < bigsi.length/1000; i++) {

        const dir = zeroPad(i, 3)

        try {
            fs.mkdirSync(`bigsiBlobs/${dir}`, { recursive: true } );
        } catch (e) {
            console.log('Cannot create folder ', e);
        }
    }

    for (let i = 0; i < bigsi.length; i++) {

        const bits = bitStringToArrayBuffer(bigsi[i])

        const paddedRowNum = zeroPad(i, 6)
        const rowPath = `${paddedRowNum.slice(0,3)}/${paddedRowNum.slice(3,6)}`
        const path = `bigsiBlobs/${rowPath}.bin`


        fs.writeFileSync(path, Buffer.from(bits), 'binary', function (err) {
            if(err) {
                console.log(err)
            } else {
                console.log(`BIGSI written to: ${path}`);
                }
        });
    }
}

function readBigsiRowArrBuf(){

    var queryRows = [];
    while(queryRows.length < 100){
        var r = Math.floor(Math.random() * 230) + 1;
        if(queryRows.indexOf(r) === -1) queryRows.push(r);
    }
    console.log('queryRows: ', queryRows);

    let bucketHits = (new BitSet).flip()
    //console.log('init bs:', bucketHits.toString())
    for (let i = 0; i < queryRows.length; i++){
        const paddedRowNum = zeroPad(queryRows[i], 6)
        const rowPath = `${paddedRowNum.slice(0,3)}/${paddedRowNum.slice(3,6)}`
        const inputFile = `bigsiBlobs/${rowPath}.bin`

        let bigsiRowBuf = fs.readFileSync(inputFile)
        console.log(bigsiRowBuf)
        const bs = new BitSet(bigsiRowBuf)
        console.log(bs.toString(16))
        bucketHits = bucketHits.and(bs)
        //console.log(`row: ${queryRows[i]}`, bs.toString(), '\n\n',  bucketHits.toString())
    }

    const bucket = bucketHits.toArray()
    console.log(bucket)
}

async function makeBinaryBigsi(){
    bigsi = []
    for (let i = 0; i < 300000; i++){
        const bs = BitSet.Random(230);
        bigsi.push(bs.toString(16))
    }

    const bigsiJSON = JSON.stringify(bigsi, null, 4);

    const outputFilename = 'binaryBigsitest.json'
    await fs.writeFile(outputFilename, bigsiJSON, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI written to: ${outputFilename}`);
            }
    });
}

function queryBigsi(){
    const binaryBigsi = require('./binaryBigsitest.json')

    var queryRows = [];
    while(queryRows.length < 100){
        var r = Math.floor(Math.random() * 230) + 1;
        if(queryRows.indexOf(r) === -1) queryRows.push(r);
    }
    console.log('queryRows: ', queryRows);

    let bucketHits = (new BitSet).flip()
    //console.log('init bs:', bucketHits.toString())
    for (let i = 0; i < queryRows.length; i++){
        const bs = BitSet.fromHexString(binaryBigsi[queryRows[i]])
        bucketHits = bucketHits.and(bs)
        //console.log(`row: ${queryRows[i]}`, bs.toString(), '\n\n',  bucketHits.toString())
    }

    const bucket = bucketHits.toArray()
    console.log(bucket)
}

function bitsetTest(){
    bstring1 = "1110"
    bstring2 = "10111"

    bs1 = new BitSet(bstring1)
    bs2 = new BitSet(bstring2)

    const hits = bs1.and(bs2).toArray()
    const hitsBucketNums = hits.map( value => 4 - value )

    console.log(hitsBucketNums)

    console.log(bs1.toString())
    console.log(bs2.toString())
    console.log(bs1.and(bs2).toString())
    console.log(bs1.and(bs2).toArray())
}

function bucketSizeTest(){
    const numBuckets = 10
    const overhang = 300_000
    const seqSize = 275_348_123
    const bucketSize = Math.round((seqSize + (numBuckets - 1 ))/numBuckets)
    
    for (let bucketNum = 0; bucketNum < numBuckets; bucketNum++){
        const bucketStart = Math.max((bucketSize-overhang)*bucketNum, 0)
        const bucketEnd = Math.min(bucketStart + bucketSize, seqSize)
        const bucketLength = bucketEnd - bucketStart

        console.log(bucketNum, bucketStart, bucketEnd, bucketLength)
    }
}

async function compareSplitToHexBigsi(){
    // pick some random rows

    const queryRows = [];
    while(queryRows.length < 10){
        const r = Math.floor(Math.random() * 299999);
        if(queryRows.indexOf(r) === -1) queryRows.push(r);
    }

    // import hexbigsi and some row files
    const hexBigsi = require('./data/hg38_chr1and2_hexbigsi.json')
    const root = './data/hg38_chr1and2_binary/'
    // extract hexbigsi row and row file
    for (let i = 0; i < queryRows.length; i++){
        const hexRow = BitSet.fromHexString(hexBigsi[queryRows[i]]).toString()

        const rowPath = commonFunc.makeBigsiRowPath(queryRows[i], root)

        let bigsiRowBuf = await fs.promises.readFile(rowPath, (err, data) => {
                if (err) {
                    console.error(err)
                    return
                }
            })

        const array = new Uint8Array(bigsiRowBuf);
        const bs = new BitSet(array)
        const binaryRow = bs.toString()
        
        console.log(queryRows[i], hexRow, binaryRow)
    }
}

async function testMakeBinaryBigsi(){
    const bitstr = '1111011'
    const bits = commonFunc.bitStringToArrayBuffer(bitstr)
    const path = './data/testMakeBinaryBigsi.bin'

    fs.writeFileSync(path, Buffer.from(bits), 'binary', function (err) {
        if(err) {
            console.log(err)
        } else {
            console.log(`BIGSI row ${row} written to: ${path}`);
            }
    });


    let bigsiRowBuf = await fs.promises.readFile(path, (err, data) => {
            if (err) {
                console.error(err)
                return
            }
        })

    const array = new Uint8Array(bigsiRowBuf);
    const bs = new BitSet(array)
    const binaryRow = bs.toString()

    console.log(bitstr, binaryRow)
}

async function testTwoBinaryBigsi(root){

    const queryRows = [91642, 233548, 58569];

    for (let i = 0; i < queryRows.length; i++){

        const rowPath = commonFunc.makeBigsiRowPath(queryRows[i], root)

        let bigsiRowBuf = await fs.promises.readFile(rowPath, (err, data) => {
                if (err) {
                    console.error(err)
                    return
                }
            })

        const array = new Uint8Array(bigsiRowBuf);
        const bs = new BitSet(array)
        const binaryRow = bs.toString()
        
        console.log(queryRows[i], binaryRow)
    }
}

async function main(){
    //const A = [[0, 1], [1,0], [0,1]]
    //const matrixA = matrix(A)
    //const transA = matrixA.trans()
    //console.log(matrixA(), transA)
}


main()


