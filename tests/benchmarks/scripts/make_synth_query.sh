BIGSI_SIZES=( 1200 )
#BIGSI_SIZES=( 100 400 800 1200 1600 2000 2500 3000)
#BIGSI_SIZES=( 200 400 800 1600 2000 2500 3000)
QUERY_ENDS=( 1000 2000 4000 8000 16000 32000 64000 128000 256000 300000 )

for query_end in ${QUERY_ENDS[@]};
do
	# handle 100M bigsi
	samtools faidx seqs/synthetic/performance/synth_100.fasta synthetic_seq:1-$query_end > seqs/synthetic/performance/synth_query_100_$query_end.fasta
	samtools faidx seqs/synthetic/performance/synth_query_100_$query_end.fasta
	for size in ${BIGSI_SIZES[@]};
	do
		samtools faidx seqs/synthetic/performance/synth_$size.fasta synthetic_seq_0:1-$query_end > seqs/synthetic/performance/synth_query_${size}_${query_end}.fasta
		samtools faidx seqs/synthetic/performance/synth_query_${size}_${query_end}.fasta
	done
done
