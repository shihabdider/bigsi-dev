#BIGSI_SIZES=( 100 200 400 800 1600 2000 2500 3000)
BIGSI_SIZES=( 1200 )

for size in ${BIGSI_SIZES[@]};
do
    time ../../../bin/query_bigsi.js \
        -r ../../../bigsis/synth_$size.bigsi.bin \
        -q seqs/synthetic/performance/synth_query_$size.fasta \
        -c ../../../bigsi.query.synth_$size.config.json \
        > metrics/synth_bigsi_size.txt
done
