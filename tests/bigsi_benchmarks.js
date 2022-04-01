const ncbi = require('bionode-ncbi')

ncbi.fetch('assembly', 'hg38').on('data', console.log)

