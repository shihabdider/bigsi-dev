const makeBigsi = require('./make_bigsi.js')
const commonFunc = require('./common_func.js')

async function main(){
    const { argv } = require('yargs')
        .scriptName('make multi-genome bigsi')
        .usage('Usage: $0 --dir (path to ref fastas and fais) --output (output path)')
        .option('dir', {
            alias: 'd',
            describe: 'Directory where reference sequences and their indices are stored (in fasta/fai format)',
            demandOption: 'Valid directory is required',
            type: 'string',
            nargs: 1,
        })
        .option('output', {
            alias: 'o',
            describe: 'Binary output file of BIGSI',
            default: 'bigsis/output_bigsi.json',
            type: 'string',
            nargs: 1,
        })

        console.log('Getting sequences in', argv.dir)
        const genomeSeqs = await makeBigsi.getMultiGenomeSeqs(argv.dir)
        console.log(`Finished extracting genomic sequences, building bigsi...`)
        const bigsi = await makeBigsi.buildBigsi(genomeSeqs)
        console.log(`Bigsi built!`)

        const memoryUsed = process.memoryUsage().heapUsed / 1024 / 1024;
        console.log(`Process uses ${memoryUsed}`)

        const u16IntRows = makeBigsi.bigsiToInts(bigsi, 16)
        const binaryDumpBigsi = makeBigsi.makeBinaryDumpBigsi(u16IntRows)
        console.log(`Converted bigsi matrix to binary TypedArray format, writing to file...`)
        makeBigsi.writeBinaryDumpBigsi(binaryDumpBigsi, `${argv.output}_bdump.bin`)
            
        const bucketToGenome = await makeBigsi.makeBucketToGenome(argv.dir)
        makeBigsi.writeBucketMapToJSON(bucketToGenome, `${argv.output}_bucket_map.json`)
}


main()
