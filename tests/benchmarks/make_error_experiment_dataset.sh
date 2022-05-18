#for i in {1..20}; 
#do
#	python make_benchmark_seqs.py -i hg38_refseq_acn.txt -o seqs/hg38.experiment.$i.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai -n 4;
#done
python make_benchmark_seqs.py -i hg38_refseq_acn.txt -o seqs/hg38.random.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai -n 40;
