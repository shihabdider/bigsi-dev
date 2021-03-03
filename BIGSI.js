const { IndexedFasta, BgzipIndexedFasta } = require('@gmod/indexedfasta')
const kmer = require('kmer.js')
const murmur = require('murmurhash-js')
const { BloomFilter } = require('bloom-filters')
const matrix = require('matrix-js')
const fs = require('fs')
const cdf = require('binomial-cdf');
 

// => 0.699

//const chr1 = new IndexedFasta({
//    path: 'Homo_sapiens.GRCh38.dna.chromosome.1.fa',
//    faiPath: 'Homo_sapiens.GRCh38.dna.chromosome.1.fa.fai',
//    chunkSizeLimit: 50000000
//});

async function testIndexedFasta(){

    // get the first 10 bases of a sequence from the file.
    // coordinates are UCSC standard 0-based half-open
    const chr1Region = await chr1.getSequence('1', 0, 10);
    // chr1Region is now a string of bases, 'ACTG...'

    // get a whole sequence from the file
    //const chr1Bases = await chr1.getSequence('1');

    // get object with all seq lengths as { seqName => length, ... }
    const allSequenceSizes = await chr1.getSequenceSizes();

    // get the size of a single sequence
    const chr1Size = await chr1.getSequenceSize('1');

    // get an array of all sequence names in the file
    const seqNames = await chr1.getSequenceNames();

    console.log(chr1Region, allSequenceSizes, chr1Size, seqNames)
}

async function splitIntoBuckets(){
    const chr1Size = await chr1.getSequenceSize('1');
    const numBuckets = 10
    const bucketSize = Math.floor(chr1Size/numBuckets)
    const overhang = 300000
    var buckets = []
    for (i = 1; i <= numBuckets; i++){
        //Handle exception for 1st bucket
        if (i == 1){
            const bucket1Region = await chr1.getSequence('1', 0, bucketSize);
            buckets.push(bucket1Region)
            //console.log(0, bucketSize, bucket1Region.length)
        } else {
            const bucketStart = ((i-1)*(bucketSize-overhang))+1
            const bucketEnd = Math.min((bucketStart + bucketSize), chr1Size)
            const bucketRegion = await chr1.getSequence('1', bucketStart, bucketEnd);
            buckets.push(bucketRegion)
            //console.log(bucketStart, bucketEnd, bucketRegion.length)
        }
    }
    return buckets
}

class Minimizer {
    constructor(hash, window_pos) {
        this.hash = hash
        this.window_pos = window_pos
    }

    isLessThan(Minimizer) {
        if (this.hash < Minimizer.hash) {
            return 1
        } else if (this.hash == Minimizer.hash) {
            if (this.window_pos < Minimizer.window_pos) {
                return 1
            } else { return 0 }
        } else { return 0 }
    }

    isEqual(Minimizer) {
       if (this.hash == Minimizer.hash && this.window_pos == Minimizer.window_pos) {
           return 1
       } else { return 0 }
    }

    isNotEqual(Minimizer) {
        if (this.hash != Minimizer.hash || this.window_pos != Minimizer.window_pos) {
            return 1
        } else { return 0 }
    }
}

//const a = new Minimizer(1, 1)
//const b = new Minimizer(2, 2)
//const c = new Minimizer(1, 1)

//console.log('Is Equal', a.isEqual(b), a.hash, b.hash)
//console.log('Less than', a.isLessThan(c), a.hash, c.hash)
//console.log('Not equal', a.isNotEqual(a), a.hash)

//const hash = murmur.murmur3('ACTGACTGACTGACTG')
//const hash2 = murmur.murmur3('ACTGACTGACTGACTG')
//console.log(hash, hash2)

function bucketToMinimizers(bucketSequence){

    const kmerSize = 16
    const windowSize = 100
    const seed = 42

    let bucketMinimizers = []
    let bucketMinimizersDeque = [] // array of {hash, offset}
    for (i = 0; i < (bucketSequence.length - kmerSize + 1); i++){
        let currentWindowIndex = i - windowSize + 1
        let kmerHash = murmur.murmur3(bucketSequence.slice(i,i+kmerSize), seed)
        //console.log(kmerHash, bucketSequence.slice(i,i+kmerSize))

        while (bucketMinimizersDeque.length != 0 && bucketMinimizersDeque[0].offset <= i - windowSize){
            bucketMinimizersDeque.shift()
        }

        while (bucketMinimizersDeque.length != 0 && bucketMinimizersDeque.slice(-1)[0].hash >= kmerHash)
        {
            bucketMinimizersDeque.pop()
        }

        bucketMinimizersDeque.push({'hash':kmerHash, 'offset':i})
        //console.log('deque', bucketMinimizersDeque.slice(-1)[0])

        if (currentWindowIndex >= 0){
            if ( bucketMinimizers.length == 0 || bucketMinimizers.slice(-1)[0] != bucketMinimizersDeque[0].hash )
            {
                bucketMinimizers.push(bucketMinimizersDeque[0].hash)
                //console.log('bucketMinimizer', bucketMinimizers.slice(-1)[0])
            }
        }
    }
    return bucketMinimizers
}

function getQueryIndices(queryBloomFilter){
    return queryBloomFilter.reduce((indices, number, index) => {
        if (number != 0) indices.push(index)
        return indices
    }, [])

}

async function buildBIGSI(){
    bucketSequences = await splitIntoBuckets()
    console.log('num buckets', bucketSequences.length)

    let bucketsBloomFilters = []
    for (j = 0; j < bucketSequences.length; j++){
    //for (j = 0; j < 2; j++){
        let bucketSequence = bucketSequences[j]
        //console.log(bucketSequence.length)
        let bucketMinimizers = bucketToMinimizers(bucketSequence)
        console.log('minimizers per bucket', bucketMinimizers.length)

        let bucketBloomFilter = new BloomFilter(300000, 1)
        for (i=0; i < bucketMinimizers.length; i++){
            bucketBloomFilter.add(bucketMinimizers[i].toString())
        }
        console.log('bloom filter', bucketBloomFilter.length)
        //console.log('bloom filter error', bucketBloomFilter.rate())
        //console.log('bloom filter check', bucketBloomFilter.has(106131787))
        //console.log('bloom filter check', bucketBloomFilter.has('test_fail'))
        bucketsBloomFilters.push(bucketBloomFilter)
        //console.log( pos(bucketBloomFilter._filter) )
        //console.log( bucketBloomFilter._filter.reduce((a, b) => a + b, 0) )

    }
    console.log('all bloom filters', bucketsBloomFilters.length)
    var M = matrix(bucketsBloomFilters.map(bloomFilterObj => bloomFilterObj._filter))
    let BIGSI = matrix(M.trans())

    var filename = 'bigsi.txt';
    var str = JSON.stringify(BIGSI(), null, 4);

    fs.writeFile(filename, str, function(err){
        if(err) {
            console.log(err)
        } else {
            console.log('File written!');
        }
    });

    //return BIGSI
    //console.log(BIGSI([0,5]))
    //console.log('equality', bucketsBloomFilters[0].equals(bucketsBloomFilters[3]))

}

async function buildQuery(query_path, query_faiPath){
    // Returns a list of minimizer arrays for each fragment
    const querySeqObj = new IndexedFasta({
        path: query_path,
        faiPath: query_faiPath,
        chunkSizeLimit: 50000000
    });

    const fragment_size = 2500
    let queryAllFragmentsMinimizers = []
    
    const query_size = await querySeqObj.getSequenceSize('1')

    if (query_size <= 300000 && query_size >= 5000){

        for (n=0; n < Math.floor(query_size/fragment_size); n++){
            let start = n*fragment_size
            let end = Math.min(start + fragment_size, query_size)
            //const [start, end] = [115285917, 117338249]
            const queryFragmentSeq = await querySeqObj.getSequence('1', start, end);
            const queryFragmentMinimizers = bucketToMinimizers(queryFragmentSeq)
            queryAllFragmentsMinimizers.push(queryFragmentMinimizers)
        }

        let queryAllBloomFilters = []
        for (j=0; j < queryAllFragmentsMinimizers.length; j++){
            let queryFragmentMinimizers = queryAllFragmentsMinimizers[j]

            let queryBloomFilter = new BloomFilter(300000, 1)
            for (k=0; k < queryFragmentMinimizers.length; k++){
                queryBloomFilter.add(queryFragmentMinimizers[k].toString())
            }
            queryAllBloomFilters.push(queryBloomFilter)
        }


        //console.log('query bloom filter', queryAllBloomFilters.length)
        return queryAllBloomFilters
    }

}

async function queryBIGSI(bigsi_path, query_path, query_faiPath){
    let rawdata = fs.readFileSync(bigsi_path);
    let parsed = JSON.parse(rawdata);
    let BIGSI = matrix(parsed)
    //console.log(BIGSI(1));

    const queryAllBloomFilters = await buildQuery(query_path, query_faiPath)
    console.log('number of query BFs: ', queryAllBloomFilters.length)

    let BIGSI_hits = {}

    for (numBuckets=0; numBuckets < BIGSI(0).length; numBuckets++){
        //console.log(numBuckets)
        BIGSI_hits[numBuckets] = 0
    }

    for (m = 0; m < queryAllBloomFilters.length; m++){
        const queryIndices = getQueryIndices(queryAllBloomFilters[m]._filter)

        let queryBIGSISubmatrix = []
        for (i = 0; i < queryIndices.length; i++){
            const queryBIGSIRow = BIGSI(queryIndices[i])
            queryBIGSISubmatrix.push(queryBIGSIRow)
        }

        queryBIGSISubmatrix = matrix(queryBIGSISubmatrix).trans()
        
        for (bucketId = 0; bucketId < queryBIGSISubmatrix.length; bucketId++){
            rowsum = queryBIGSISubmatrix[bucketId].reduce((a, b) => a + b, 0)
            const score = cdf(rowsum,queryBIGSISubmatrix[bucketId].length, 0.82)
            if (score >= 0.999){
                BIGSI_hits[bucketId] += 1
            }
        }
        
    }
    
    //console.log(BIGSI_hits)
    return BIGSI_hits
}

//console.log(hash, hash2)
//testIndexedFasta();
//buildBIGSI()
async function main(){
    // path: 'Homo_sapiens.GRCh38.dna.chromosome.1.fa',
    // faiPath: 'Homo_sapiens.GRCh38.dna.chromosome.1.fa.fai',
    // path: 'Homo_sapiens.GRCh38.dna.chromosome.1.split.100Kmer.007.fa',
    // faiPath: 'Homo_sapiens.GRCh38.dna.chromosome.1.split.100Kmer.007.fa.fai',
    const path = 'Homo_sapiens.GRCh38.dna.chromosome.1.split.100Kmer.007.fa'
    const faiPath = 'Homo_sapiens.GRCh38.dna.chromosome.1.split.100Kmer.007.fa.fai'
    const hits = await queryBIGSI('bigsi.txt', path, faiPath)
    console.log(hits)
}

main()
