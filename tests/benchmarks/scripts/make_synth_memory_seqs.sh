#!/bin/zsh
#BIGSI_SIZES=( 100 400 800 1200 1600 2000 2500 3000)
BIGSI_SIZES=( 1200 )

# Make the synthetic reference seqs
for size in ${BIGSI_SIZES[@]};
do
    python scripts/make_synthetic_seq.py $size \
        seqs/synthetic/performance/synth_$size.fasta;
    samtools faidx seqs/synthetic/performance/synth_$size.fasta \
	&
done
wait
