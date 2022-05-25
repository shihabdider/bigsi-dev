BIGSI_SIZES=( 100 200 400 800 1600 2000 2500 3000)

# Make the synthetic reference BIGSIs
for size in ${BIGSI_SIZES[@]};
do
    node ../../bigsi.js \
        -r seqs/synthetic/performance/synth_$size.fasta \
        -o bigsis/synth_$size.bigsi \
	&
done

