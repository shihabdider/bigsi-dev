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
