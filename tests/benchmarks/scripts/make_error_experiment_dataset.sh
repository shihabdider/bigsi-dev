for i in {1..100}; 
do
	python scripts/make_benchmark_seqs.py \
		-i scripts/hg38_refseq_acn.txt \
		-o seqs/hg38/simulation/substitution_rate/experiment_$i/hg38_sample.fasta \
		-f ~/Research/seqs/hg38.fna \
		-x ~/Research/seqs/hg38.fna.fai \
		-n 4;
done
#python make_benchmark_seqs.py -i hg38_refseq_acn.txt -o seqs/hg38.random.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai -n 40;
#python make_benchmark_seqs.py -i synth_acn.txt -o seqs/synthetic_seq_300M.random.fasta -f seqs/synthetic_seq_300M.fasta -x seqs/synthetic_seq_300M.fasta.fai -n 1000;
