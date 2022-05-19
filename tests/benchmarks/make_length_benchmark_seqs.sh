#usage: make_benchmark_seqs.py [-h] -i IDENTIFIERS [-l LENGTH] [-n NUM] -o OUTPUT
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 5000 -n 40 -o seqs/hg38_5KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 7000 -n 40 -o seqs/hg38_7KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 20000 -n 40 -o seqs/hg38_20KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 40000 -n 40 -o seqs/hg38_40KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 80000 -n 40 -o seqs/hg38_80KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 160000 -n 40 -o seqs/hg38_160KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 200000 -n 40 -o seqs/hg38_200KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 250000 -n 40 -o seqs/hg38_250KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -l 300000 -n 40 -o seqs/hg38_300KB.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai

python make_benchmark_seqs.py -i synth_acn.txt -l 1000 -n 1000 -o seqs/synthetic_seq_300M.1000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 2000 -n 1000 -o seqs/synthetic_seq_300M.2000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 3000 -n 1000 -o seqs/synthetic_seq_300M.3000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 4000 -n 1000 -o seqs/synthetic_seq_300M.4000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 5000 -n 1000 -o seqs/synthetic_seq_300M.5000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 10000 -n 1000 -o seqs/synthetic_seq_300M.10000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 20000 -n 1000 -o seqs/synthetic_seq_300M.20000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 40000 -n 1000 -o seqs/synthetic_seq_300M.40000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 80000 -n 1000 -o seqs/synthetic_seq_300M.80000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 160000 -n 1000 -o seqs/synthetic_seq_300M.160000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 200000 -n 1000 -o seqs/synthetic_seq_300M.200000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 250000 -n 1000 -o seqs/synthetic_seq_300M.250000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
python make_benchmark_seqs.py -i synth_acn.txt -l 300000 -n 1000 -o seqs/synthetic_seq_300M.300000.fasta -f ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta -x ~/Research/bigsi-dev/tests/benchmarks/seqs/synthetic_seq_300M.fasta.fai
