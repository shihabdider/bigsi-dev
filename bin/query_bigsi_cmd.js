const queryBigsi = require('./query_bigsi.js')
const utils = require('./utils.js')

function printResults(filteredBigsiHits, binMapPath) {
    // import the associated map json file as dictionary
    const binMap = require(binMapPath)
    // iterate through filtered hits keys
    for (key in filteredBigsiHits) {
        const {refName, bucketStart, bucketEnd}  = binMap[key]
        const outputString = `${refName} ${bucketStart} ${bucketEnd}`
        console.log(outputString) 
    }
}

async function run() {
    if (require.main === module) {
        const { argv } = require('yargs')
            .scriptName('query_bigsi')
            .usage('Usage: $0 --query (path to query fasta) OR --querySeq --bigsi (path to BIGSI file) --config (path to bigsi config)')
            .option('querySeq', {
                alias: 's',
                describe: 'Query sequence',
                type: 'string',
                nargs: 1,
            })
            .option('query', {
                alias: 'q',
                describe: 'Path to fasta file of query sequence',
                type: 'string',
                nargs: 1,
            })
            .option('bigsi', {
                alias: 'b',
                describe: 'Path to BIGSI file',
                demandOption: 'BIGSI file is required',
                type: 'string',
                nargs: 1,
            })
            .option('subrate', {
                alias: 'e',
                describe: 'Substitution rate threshold',
                demandOption: 'Substitution rate threshold is required',
                type: 'number',
                nargs: 1,
            })
            .option('config', {
                alias: 'c',
                describe: 'Path to BIGSI config file',
                demandOption: 'BIGSI config file is required',
                type: 'string',
                nargs: 1,
            })

        const bigsiPath = argv.bigsi
        const bigsiConfigPath = argv.config
        const binMapPath = bigsiPath.slice(0, -4) + '_bucket_map.json'
        if (argv.querySeq) { // passing query as string
            const hits = await queryBigsi.main(argv.querySeq, bigsiPath, bigsiConfigPath, argv.subrate);
            printResults(hits, binMapPath)
        } else { // ...or else as a file
            const fai = `${argv.query}.fai`
            const query = await utils.loadFasta(argv.query, fai)
            const seqNames = await query.getSequenceList()
            //console.log('Query sequence name: ', seqNames)
            for (seqName of seqNames) {
                const querySeq = await query.getSequence(seqName)
                const hits = await queryBigsi.main(querySeq, bigsiPath, bigsiConfigPath);
                printResults(hits, binMapPath)
            }
        }
    }
}

run()

