#BIGSI_SIZES=( 100 400 800 1200 1600 2000 2500 3000)
#SIZES_MB=( 0.05 0.1 0.2 0.4 0.8 1.6 3 )
SIZES_MB=( 16 )
query_sizes=( 20000 40000 80000 160000 )

function make_synth_seq() {
    # Make the synthetic reference seqs
    local output_dir=$1
    for size in ${SIZES_MB[@]};
    do
        python3 scripts/make_synthetic_seq.py ${size} \
            seqs/synthetic/${1}/${size}.fasta;
        samtools faidx seqs/synthetic/${1}/${size}.fasta \
        &
    done
    wait
}

function make_synth_query() {
    local num_queries_per_seq=1000
    local dir=$1
    local ref_size=16

    for query_size in ${query_sizes[@]};
    do
        local output=seqs/synthetic/${1}/${ref_size}_${query_size}_query_000.fasta
        python3 scripts/make_benchmark_seqs.py \
            -i scripts/synth_acn.txt \
            -l ${query_size} \
            -n ${num_queries_per_seq} \
            -o ${output}\
            -f seqs/synthetic/$1/${size}.fasta \
            -x seqs/synthetic/$1/${size}.fasta.fai;
        samtools faidx ${output};
    done
}

function mutate_synth_seq() {
    sub_rates=( 001 002 003 004 005 006 007 008 009 010 )
    local ref=$1
    local N=12
    local ref_size=16
    for query_size in ${query_sizes[@]};
    do
        for i in ${!sub_rates[@]};
        do
            ((job=job%N)); ((job++==0)) && wait
            output=seqs/synthetic/${1}/${ref_size}_${query_size}_query_${sub_rates[i]}.fasta;
            python3 scripts/mutate_seqs.py \
                -i seqs/synthetic/${1}/${ref_size}_${query_size}_query_000.fasta \
                -r ${sub_rates[i]} \
                -o ${output};
            samtools faidx ${output} \
            &
        done
    done
}

make_synth_seq jaccard;
make_synth_query jaccard;
mutate_synth_seq jaccard;
