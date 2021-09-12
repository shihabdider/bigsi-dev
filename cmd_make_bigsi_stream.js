const makeBigsi = require('./make_bigsi.js')
const commonFunc = require('./common_func.js')

async function main(){
    const { argv } = require('yargs')
        .scriptName('make bigsi stream')
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
        .option('output', {
            alias: 'o',
            describe: 'Output file (JSON) of BIGSI',
            default: 'bigsis/output_bigsi.json',
            type: 'string',
            nargs: 1,
        });

    const seq = await commonFunc.loadFasta(argv.fasta, argv.fai)
    const seqName = 2*10**7    // 
    
    const bigsis = await makeBigsi.makeGenomeBigsis(seq, seqSizeThreshold)
    console.log('num bigsis', bigsis.length)
    const bigsi = await makeBigsi.mergeBigsis(bigsis)

    const hexBigsi = makeBigsi.makeHexBigsi(bigsi)
    makeBigsi.writeBigsiToJSON(hexBigsi, `${argv.output}_hex.json`)

    //const bitBigsi = makeBigsi.makeBitBigsi(bigsi)
    //makeBigsi.writeBigsiToJSON(bitBigsi, `${argv.output}_bit.json`)

    //const bitBigsi = makeBigsi.makeBitBigsi(bigsi)
    //const numDirs = 300
    //makeBigsi.writeBigsiArrayBuf(numDirs, argv.dir, bitBigsi)

    //const bucketToPosition = await makeBigsi.makeBucketToPosition(seq, seqSizeThreshold)
    //makeBigsi.writeBucketMapToJSON(bucketToPosition, `${argv.output}_bucket_map.json`)

}

main()
