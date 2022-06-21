query_size=300000
function jaccard_benchmark() {
    local sub_rates=( 000 001 002 003 004 005 006 007 008 009 010 )
    local ref_sizes=( 16 )
    window_size=$1
    N=12
    # iterate through all refs
    for size in ${ref_sizes[@]};
    do
        # iterate through all mutated queries for that ref
        for rate in ${sub_rates[@]};
        do
            ((b=b%N)); ((b++==0)) && wait
            ref="seqs/synthetic/jaccard/${size}.fasta"
            query="seqs/synthetic/jaccard/${size}_${query_size}_query_${rate}.fasta"
            node ~/Research/bigsi-dev/tests/jaccard_test.js \
                -r ${ref} -q ${query} -w ${window_size} > \
                metrics/jaccard/${size}_${query_size}_${rate}_w${window_size}.txt \
        &
        done
        echo "Benchmarks for $size reference complete!"
    done
}

jaccard_benchmark 25;
jaccard_benchmark 50;
jaccard_benchmark 100;
