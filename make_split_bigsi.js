/**
 * @param { matrix } bigsi - bigsi matrix
 *
 * @returns { Array of binary strings } bitBigsi
 */
function makeBitBigsi(bigsi){
    let bitBigsi = []
    const numRows = bigsi.size()[0]

    for (let i = 0; i < numRows; i++){
        const ithRow = bigsi(i)
        const binaryString = ithRow.join('')
        bitBigsi.push(binaryString)
    }

    return bitBigsi
}
/**
 * @param { number } numDirs - number of parent directories for storing bigsi 
 * rows
 * @param { string } root - parent directory for storing bigsi dirs
 * rows
 *
 * @returns - directories for storing bigsi rows, numRows/directory
 */
function makeBigsiDirs(numDirs, root='bigsiBinaryFiles'){

    for (let i = 0; i < numDirs; i++) {

        const dir = commonFunc.zeroPad(i, 3)

        try {
            fs.mkdirSync(`${root}/${dir}`, { recursive: true } );
        } catch (e) {
            console.log('Cannot create folder ', e);
        }
    }
}
/**
 * @param { array of binary strings } bitBigsi - bigsi matrix as an array of binary 
 * strings (string/row)
 *
 * @returns - one binary file per row/element of the bitBigsi
 */
function writeBigsiArrayBuf(numDirs, root, bitBigsi){

    makeBigsiDirs(numDirs, root)

    const numRows = bitBigsi.length
    for (let row = 0; row < numRows; row++) {

        const path = commonFunc.makeBigsiRowPath(row, root)
        const bits = commonFunc.bitStringToArrayBuffer(bitBigsi[row])

        fs.writeFileSync(path, Buffer.from(bits), 'binary', function (err) {
            if(err) {
                console.log(err)
            } else {
                console.log(`BIGSI row ${row} written to: ${path}`);
                }
        });
    }
}

