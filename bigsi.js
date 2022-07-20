/* Command line tool for building bigsi matrix of a fasta sequence
 */

const makeBigsi = require('./bin/make_bigsi.js')
const writeBigsi = require('./bin/write_bigsi.js')
const utils = require('./bin/utils.js')

async function main(){
    const { argv } = require('yargs')
        .scriptName('bigsi')
        .usage('Usage: $0 --ref (path to ref fasta) --output (output path)')
        .option('ref', {
            alias: 'r',
            describe: 'Fasta file of reference sequence',
            demandOption: 'Fasta file is required',
            type: 'string',
            nargs: 1,
        })
        .option('output', {
            alias: 'o',
            describe: 'Output file of BIGSI',
            default: 'bigsis/output',
            type: 'string',
            nargs: 1,
        })

    const fai = `${argv.ref}.fai`
    const fasta = await utils.loadFasta(argv.ref, fai)
    console.log('Sequence loaded...')

    //const [bigsi, bigsiDims] = await makeBigsi.main(fasta)
    //const bigsiInts = writeBigsi.bitstringsToInts(bigsi)
    //const binaryBigsi = writeBigsi.makeBinaryBigsi(bigsiInts)
    //console.log(`Converted bigsi matrix to binary TypedArray format, writing to file...`)
    //writeBigsi.writeBinaryBigsi(binaryBigsi, `${argv.output}.bin`)
    //writeBigsi.writeQueryConfigToJSON(bigsiDims, `${argv.output}_query_config.json`)

    const seqToPositionMap = await writeBigsi.makeSeqToPositionMap(fasta)
    writeBigsi.writeSeqMapToJSON(seqToPositionMap, `${argv.output}_bucket_map.json`)
}

main()
