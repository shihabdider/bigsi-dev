const { BloomFilter } = require('bloom-filters')
// create a Bloom Filter with a size of 10 and 4 hash functions
let filter = new BloomFilter(10, 1)
// insert data
filter.add('102341044')
filter.add('102341045')

// lookup for some data
console.log(filter.has('102341044')) // output: true
console.log(filter.has('102033044')) // output: false

// print the error rate
console.log(filter.rate())
console.log(filter._filter)
console.log(filter._length)

console.log( filter._filter.reduce((a, b) => a + b, 0) )
