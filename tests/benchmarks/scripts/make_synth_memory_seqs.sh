BIGSI_SIZES=( 100 200 400 800 1600 2000 2500 3000)

# Make the synthetic reference seqs
for size in ${BIGSI_SIZES[@]};
do
    python scripts/make_synthetic_seq.py $size \
        seqs/synthetic/performance/synth_$size.fasta 
done

# Make the synthetic query seqs
for size in ${BIGSI_SIZES[@]};
do
    python scripts/make_benchmark_seqs.py \
		-i scripts/synth_acn.txt \
		-o seqs/synthetic/performance/synth_query_$size.fasta \
		-f seqs/synthetic/performance/synth_$size.fasta \
		-x seqs/synthetic/performance/synth_$size.fasta.fai \
		-n 1;
done

