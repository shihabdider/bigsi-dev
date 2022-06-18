HG38_QUERY_LEN=hg38/simulation/query_length
QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
for i in ${QUERY_LENGTHS[@]};
do
	for j in {1..100};
	do
		python scripts/make_benchmark_seqs.py \
			-i scripts/hg38_refseq_acn.txt \
			-l ${i} \
			-n 4 \
			-o seqs/${HG38_QUERY_LEN}/experiment_${j}/${i}.fasta \
			-f ~/Research/seqs/hg38.fna \
			-x ~/Research/seqs/hg38.fna.fai
	done
done

