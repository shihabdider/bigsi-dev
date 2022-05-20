sizes=( 1000 2000 3000 4000 ) 
#5000 7000 10000 20000 40000 80000 160000 100000 200000 250000 300000 )
for size in ${sizes[@]}; do
	for i in {1..20}; 
	do
		python make_benchmark_seqs.py -i hg38_refseq_acn.txt -o seqs/hg38.experiment.$i.$size.fasta -f ~/Research/seqs/hg38.fna -x ~/Research/seqs/hg38.fna.fai -n 4 -l $size;
	done
done
