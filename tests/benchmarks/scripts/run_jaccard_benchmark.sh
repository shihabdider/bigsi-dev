#query_sizes=( 5000 300000 )
query_sizes=( 20000 40000 80000 160000 )

function jaccard_benchmark() {
    local sub_rates=( 000 001 002 003 004 005 006 007 008 009 010 )
    local ref_size=16
    window_size=$1
    N=12
    # iterate through all query sizes
    for query_size in ${query_sizes[@]};
    do
        # iterate through all mutated queries for that ref
        for rate in ${sub_rates[@]};
        do
            ((b=b%N)); ((b++==0)) && wait
            ref="seqs/synthetic/jaccard/${ref_size}.fasta"
            query="seqs/synthetic/jaccard/${ref_size}_${query_size}_query_${rate}.fasta"
            node ~/Research/bigsi-dev/tests/jaccard_test.js \
                -r ${ref} -q ${query} -w ${window_size} > \
                metrics/jaccard/${ref_size}_${query_size}_${rate}_w${window_size}.txt \
        &
        done
    done
}

#jaccard_benchmark 25;
jaccard_benchmark 50;
jaccard_benchmark 100;
