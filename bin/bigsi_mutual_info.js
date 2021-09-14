const BitSet = require('bitset')
const commonFunc = require('./common_func.js')
const matrix = require('matrix-js')
const fs = require('fs');

/**
 * @param { array } column - column i of bigsi matrix
 * @param { number } x - x in {0,1} to compute marginal (e.g 0)
 *
 * @returns { number } marginalProb - marginal probability of columns
 */
function computeColumnMarginalProb(column,x){
    const length = column.length
    const sum = column.reduce((a, b) => a + b, 0)

    let marginalProb

    if (x === 0){
        marginalProb = (length-sum)/(length)
    } else {
        marginalProb = (sum)/(length)
    }

    return marginalProb
}


/**
 * @param { array of arrays } columns - columns of bigsi matrix (i.e 
 * bigsi.trans())
 *
 * @returns { array } entropies - entropies of columns
 */
function computeColumnEntropies(columns){
    const X = [0,1]

    const entropies = []
    for (column of columns){
        entropy = 0
        for (x of X){
            marginalProb = computeColumnMarginalProb(column, x)
            entropy += -1*marginalProb*Math.log2(marginalProb)
        }
        entropies.push(entropy)
    }

    return entropies
}

/**
 * @param { array of arrays } columns - column i,j of bigsi matrix
 * @param { string } xy - (x,y) in {0,1} to compute joint (e.g "01")
 *
 * @returns { number } jointProb - joint probability of columns
 */
function computeColumnJointProb(columns, xy){
    const length = columns[0].length
    let numMatches = 1

    for (let i=0; i < length; i++){
        const row = [columns[0][i], columns[1][i]].join('')
        if (row === xy){
            numMatches++
        }
    }

    const jointProb = numMatches/(length + 4)

    return jointProb
}

/**
 * @param { array of arrays } columns - column i,j of bigsi matrix
 * 
 * @returns { number } mi - mutual information between columns i and j
 */
function computeMutualInfo(columns){
    // iterate over all possible xy
    const xys = ['00', '01', '10', '11']

    let mi = 0
    for (let xy of xys){
    // get the marginal and joint probs
        const jointProb = computeColumnJointProb(columns, xy)
        const logJointProb = Math.log2(jointProb)

        const x = parseInt(xy.slice(0,1))
        const y = parseInt(xy.slice(1,2))
        const ithMarginalProb = computeColumnMarginalProb(columns[0], x)
        const ithLogMarginalProb = Math.log2(ithMarginalProb)

        const jthMarginalProb = computeColumnMarginalProb(columns[1], y)
        const jthLogMarginalProb = Math.log2(jthMarginalProb)


        const mutualInfoElem = jointProb*(logJointProb - (ithLogMarginalProb + jthLogMarginalProb))
        mi += mutualInfoElem
        //console.log(xy, jointProb, ithMarginalProb, jthMarginalProb, mutualInfoElem)
    }

    return mi
}

/** Converts hexbigsi into an matrix
 * @param { JSON } bigsi - JSON of bigsi file
 * @param { number } numBuckets - number of columns/buckets in the bigsi
 *
 * @returns { array } bigsiMatrix - bigsi as an matrix
 */ 
function makeBigsiMatrix(bigsi, numBuckets){
    const bigsiArr = []
    const rowLength = numBuckets
    for (const row of bigsi){
        const rowArrStr = commonFunc.zeroPadBitstring((new BitSet(`0x${row}`)).toString(), rowLength).split('')
        const rowArrNum = rowArrStr.map(Number)
        bigsiArr.push(rowArrNum)
    }

    console.log('bigsi array loaded...')

    const bigsiMatrix = matrix(bigsiArr)

    return bigsiMatrix
}


/**
 * @param { JSON } bigsi - JSON of bigsi file
 * @param { number } numBuckets - number of columns/buckets in the bigsi
 *
 * @returns { array } bigsiMI - mutual info between all column pairs of bigsi
 */
function computeBigsiMI(bigsi, numBuckets){
    const bigsiMatrix = makeBigsiMatrix(bigsi, numBuckets)
    const numColumns = (bigsiMatrix.size())[1]
    
    const bigsiMI = []
    for (let i=0; i < numColumns/2; i++){
        for (let j=i+1; j < numColumns; j++){
            const ithColumn = matrix(bigsiMatrix([], i)).trans()
            const jthColumn = matrix(bigsiMatrix([], j)).trans()

            const ijColumns = [ithColumn[0], jthColumn[0]]
            const mi = computeMutualInfo(ijColumns)
            bigsiMI.push([mi])
            if (mi > 0.1){
                console.log(`Tested columns ${i}, ${j}: ${mi}`)
            }
        }
    }
    
    return bigsiMI
}


/**
 * @param { JSON } bigsi - JSON of bigsi file
 * @param { number } numBuckets - number of columns/buckets in the bigsi
 *
 * @returns { array } bigsiDist - mutual info based distances between all 
 * column pairs of bigsi
 */
function computeBigsiDist(bigsi, numBuckets){
    const bigsiMatrix = makeBigsiMatrix(bigsi, numBuckets)
    const numColumns = (bigsiMatrix.size())[1]
    
    const bigsiDist = []
    for (let i=0; i < numColumns/2; i++){
        for (let j=i+1; j < numColumns; j++){
            const ithColumn = matrix(bigsiMatrix([], i)).trans()
            const jthColumn = matrix(bigsiMatrix([], j)).trans()

            const ijColumns = [ithColumn[0], jthColumn[0]]
            const mi = computeMutualInfo(ijColumns)
            const dist = mi
            bigsiDist.push([dist])
        }
    }

    const maxMI = Math.max(...bigsiDist.flat())
    for (dist in bigsiDist){
        dist[0] = maxMI - dist[0]
    }
    
    return bigsiDist
}

function computeCompression(bigsi){
    const [numRows, numColumns] = bigsi.size()
    const S1 = numRows*numColumns
    console.log('S1', S1)
    
    const columns = bigsi.trans()
    const entropies = computeColumnEntropies(columns)
    const entropiesSum = entropies.reduce((a,b) => a+b, 0)
    const S2 = numRows * entropiesSum
    console.log('S2', S2)

    // 0.724 - chr1 vs chr1
    // 0.00700 - chr1 vs chr2
    const S3 = S2 - (numRows*0.007)
    console.log('S3', S3)

    const compression_rate = S3/S1
    console.log('compression size: ', compression_rate)
}

function writeToJSON(filename, jsonContent){
    fs.writeFile(filename, jsonContent, 'utf8', function (err) {
        if (err) {
            return console.log(err);
        }

        console.log("The file was saved!");
    });
}

function writeBigsiDist(hexBigsi, filename){

    const numBuckets = 20
    const bigsiDist = computeBigsiDist(hexBigsi, numBuckets)

    const jsonContent = JSON.stringify(bigsiDist);
    writeToJSON(filename, jsonContent)
}

function main(){
    // require hexbigsi file
    const hexBigsi = require('./tests/data/c_elegans_vs_c_briggsae_hex.json')
    console.log(hexBigsi)
    //const bigsiMatrix = makeBigsiMatrix(hexBigsi, 20)
    //computeCompression(bigsiMatrix)
    writeBigsiDist(hexBigsi, './tests/data/c_elegans_c_briggsae_dist.json')
}

main()
