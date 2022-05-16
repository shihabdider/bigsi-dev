#usage: make_benchmark_seqs.py [-h] -i IDENTIFIERS [-l LENGTH] [-n NUM] -o OUTPUT
python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 5000 -n 40 -o seqs/hg38_5KB.fasta -f ~/Research/seqs/hg38.fna
python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 100000 -n 40 -o seqs/hg38_100KB.fasta -f ~/Research/seqs/hg38.fna
python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 300000 -n 40 -o seqs/hg38_300KB.fasta -f ~/Research/seqs/hg38.fna
