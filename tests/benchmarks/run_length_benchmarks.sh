#time python benchmark.py -q seqs/hg38_5KB.fasta -c hg38.office.config.json -o output/hg38_5KB -i 100
#time python benchmark.py -q seqs/hg38_7KB.fasta -c hg38.office.config.json -o output/hg38_7KB -i 100
#time python benchmark.py -q seqs/hg38_10KB.fasta -c hg38.office.config.json -o output/hg38_10KB -i 100 -l 5000
time python benchmark.py -q seqs/hg38_20KB.fasta -c hg38.office.config.json -o output/hg38_20KB -i 100 -l 10000
#time python benchmark.py -q seqs/hg38_40KB.fasta -c hg38.office.config.json -o output/hg38_40KB -i 100 -l 20000
#time python benchmark.py -q seqs/hg38_80KB.fasta -c hg38.office.config.json -o output/hg38_80KB -i 100 -l 40000
#time python benchmark.py -q seqs/hg38_160KB.fasta -c hg38.office.config.json -o output/hg38_160KB -i 100 -l 80000 
#time python benchmark.py -q seqs/hg38_200KB.fasta -c hg38.office.config.json -o output/hg38_200KB -i 100 -l 10000 
#time python benchmark.py -q seqs/hg38_250KB.fasta -c hg38.office.config.json -o output/hg38_250KB -i 100 -l 125000 
#time python benchmark.py -q seqs/hg38_300KB.fasta -c hg38.office.config.json -o output/hg38_300KB -i 100 -l 150000
