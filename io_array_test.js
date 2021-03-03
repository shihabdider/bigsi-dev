var fs = require('fs');
var arr = [[0,1], [1,0], [0,1]];
var filename = 'output.txt';
var str = JSON.stringify(arr, null, 4);

fs.writeFile(filename, str, function(err){
    if(err) {
        console.log(err)
    } else {
        console.log('File written!');
        let rawdata = fs.readFileSync('output.txt');
        let student = JSON.parse(rawdata);
        console.log(arr, student);
    }
});

