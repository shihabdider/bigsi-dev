const { BloomFilter } = require('bloom-filters')

function makeMinimizersBloomFilter(minimizers, bloomFilterSize, numHashes) {
    // adjust filter size based on number of inserted elements and desired false pos 
    // rate
    const minimizersBloomFilter = new BloomFilter(bloomFilterSize, numHashes)
    for (const minimizer of minimizers){
        minimizersBloomFilter.add(minimizer.toString())
    }
    return minimizersBloomFilter
}

function getBFlength(numBitsSet, numHashes, bfSize) {
    const length = (-1*bfSize/numHashes)*Math.log(1 - numBitsSet/bfSize)
    console.log('raw length', length)
    return Math.round(length)
}

function getBFIntersection(query, ref) {
    const intersection = []
    for (i=0; i< query._filter.length; i++) {
        if (query._filter[i] == 1) {
            intersection.push(ref._filter[i])
        }
    }

    return intersection
}

function getBFUnion(query, ref) {
    const union = []
    for (i=0; i< query._filter.length; i++) {
        if (query._filter[i] == 1 || ref._filter[i] == 1) {
            union.push(1)
        } else {
            union.push(0)
        }
    }

    return union
}

function computeNumHashes(bloomFilterSize, numElementsInserted) {
    const numHashes = Math.ceil((bloomFilterSize/numElementsInserted)*Math.log(2))
    return numHashes
}

function computeOptimalBFSize(numInserted, falsePosRate){
    return Math.ceil(-numInserted*Math.log(falsePosRate)/(Math.log(2)**2))
}

function main() {
    const query = [1, 2, 3, 4, 5, 6, 7]
    const ref = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14]
    const falsePosRate = 0.01
    const bfSize = computeOptimalBFSize(ref.length, falsePosRate)
    console.log('bfSize', bfSize)
    const numHashes = computeNumHashes(bfSize, ref.length)
    console.log('optimal num hashes', numHashes)
    const queryBF = makeMinimizersBloomFilter(query, bfSize, numHashes)
    const refBF = makeMinimizersBloomFilter(ref, bfSize, numHashes)
    const intersection = getBFIntersection(queryBF, refBF)
    const union = getBFUnion(queryBF, refBF)
    //console.log('intersection', intersection, intersection.length)
    //console.log('union', union, union.length)

    //const numBitsSetIntersection = intersection.reduce((sum, elem) => sum + elem, 0)
    const numBitsSetUnion = union.reduce((sum, elem) => sum + elem, 0)
    const unionLength = getBFlength(numBitsSetUnion, numHashes, union.length)
    console.log('union length', unionLength)
    const intersectionLength = queryBF.length + refBF.length - unionLength
    console.log('intersection length', intersectionLength)

    const bfSizeBin = computeOptimalBFSize(320000, 0.023)
}

main()


