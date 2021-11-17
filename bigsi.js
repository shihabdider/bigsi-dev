const makeBigsi = require('./bin/make_bigsi.js')
const helper = require('./bin/helper.js')

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
        .option('buckets', {
            alias: 'b',
            describe: 'Number of buckets per sequence (should be a multiple of 8)',
            default: 16,
            type: 'number',
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
    const fasta = await helper.loadFasta(argv.ref, fai)
    console.log('Sequence loaded...')

    const bigsi = await makeBigsi.main(fasta, argv.buckets)
    //console.log(`Converted bigsi matrix to binary TypedArray format, writing to file...`)
    //makeBigsi.writeBinaryBigsi(binaryBigsi, `${argv.output}.bin`)
            
    //const bucketToPosition = await makeBigsi.makeBucketToPositionMap(seq, seqSizeThreshold)
    //makeBigsi.writeBucketMapToJSON(bucketToPosition, `${argv.output}_bucket_map.json`)

}

main()
