
function jaccard_benchmark() {
    sub_rates=( 001 002 003 004 005 006 007 008 009 010 )
    ref_sizes=( 0.05 0.1 0.2 0.4 0.8 1.6 3 )
    window_size=25
    N=12
    # iterate through all refs
    for size in ${ref_sizes[@]};
    do
        # handle exact match
        ref="seqs/synthetic/jaccard/${size}.fasta"
        query="seqs/synthetic/jaccard/${size}_query.fasta"
        node ~/Research/bigsi-dev/tests/jaccard_test.js \
            -r ${ref} -q ${query} -w ${window_size} > \
            metrics/jaccard_test_${size}_000.txt

        # iterate through all mutated queries for that ref
        for rate in ${sub_rates[@]};
        do
            ((b=b%N)); ((b++==0)) && wait
            ref="seqs/synthetic/jaccard/${size}.fasta"
            query="seqs/synthetic/jaccard/${size}_query_${rate}.fasta"
            node ~/Research/bigsi-dev/tests/jaccard_test.js \
                -r ${ref} -q ${query} -w ${window_size} > \
                metrics/jaccard_test_${size}_${rate}.txt \
        &
        done
        echo "Benchmarks for $size reference complete!"
    done
}

jaccard_benchmark
