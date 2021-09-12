const makeBigsi = require('./make_bigsi.js')
const commonFunc = require('./common_func.js')

async function main(){
    const { argv } = require('yargs')
        .scriptName('make bigsi')
        .usage('Usage: $0 --fasta (path to ref fasta) --fai (path to ref fasta index) --output (output path)')
        .option('fasta', {
            alias: 'f',
            describe: 'Fasta file of reference sequence',
            demandOption: 'Fasta file is required',
            type: 'string',
            nargs: 1,
        })
        .option('fai', {
            describe: 'Fasta index file of reference sequence',
            demandOption: 'Index file is required',
            type: 'string',
            nargs: 1,
        })
        .option('name', {
            describe: 'Name of sequence',
            default: '1',
            type: 'string',
            nargs: 1,
        })
        .option('numBuckets', {
            describe: 'Number of buckets per sequence (should be a multiple of 8)',
            default: 16,
            type: 'number',
            nargs: 1,
        })
        .option('output', {
            alias: 'o',
            describe: 'Output file (JSON) of BIGSI',
            default: 'bigsis/output_bigsi.json',
            type: 'string',
            nargs: 1,
        })

    const seq = await commonFunc.loadFasta(argv.fasta, argv.fai)
    console.log('Sequence loaded...')

    const seqSizeThreshold = 2*10**6    // >seqSizeThreshold sequences only
    console.log(`Filtering sequences smaller than ${seqSizeThreshold}...`)
    
    const numBuckets = 16
    const bigsis = await makeBigsi.makeGenomeBigsis(seq, numBuckets, seqSizeThreshold)
    console.log(`Bigsis for ${bigsis.length} sequences created, merging...`)
    const bigsi = await makeBigsi.mergeBigsis(bigsis)
    console.log(`Bigsis merged!`)

    const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
    console.log(`Process uses ${memoryUsed}`)

    const hexBigsi = makeBigsi.makeHexBigsi(bigsi)
    console.log(`Converted bigsi matrix to hex format, writing to file...`)
    makeBigsi.writeBigsiToJSON(hexBigsi, `${argv.output}_hex.json`)
    
    const bucketToPosition = await makeBigsi.makeBucketToPosition(seq, seqSizeThreshold)
    makeBigsi.writeBucketMapToJSON(bucketToPosition, `${argv.output}_bucket_map.json`)

}

main()
