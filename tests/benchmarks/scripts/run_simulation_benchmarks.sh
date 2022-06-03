function pan_trog_benchmark() {
    pan_trog_dir=pan_trog/simulation
    QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    MIN_LENGTHS=( 500 1000 1500 2000 2500 5000 10000 20000 40000 80000 100000 125000 150000 )


    for j in {1..100};
    do
        mkdir -p outputs/${pan_trog_dir}/experiment_$j/;
    done

    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in {1..100};
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${pan_trog_dir}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${pan_trog_dir}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 95 \
                -l ${MIN_LENGTHS[i]} \
        &
        done
    done
}

function gorilla_benchmark() {
    gorilla_dir=gorilla/simulation
    QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    MIN_LENGTHS=( 500 1000 1500 2000 2500 5000 10000 20000 40000 80000 100000 125000 150000 )


    for j in {1..100};
    do
        mkdir -p outputs/${gorilla_dir}/experiment_$j/;
    done

    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in {1..100};
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${gorilla_dir}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${gorilla_dir}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 95 \
                -l ${MIN_LENGTHS[i]} \
        &
        done
    done
}

function query_length_benchmark() {
    HG38_QUERY_LEN=hg38/simulation/query_length
    QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    MIN_LENGTHS=( 500 1000 1500 2000 2500 5000 10000 20000 40000 80000 100000 125000 150000 )

    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in {1..100};
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 100 \
            -l ${MIN_LENGTHS[i]} \
        &
        done
    done
}

function error_benchmark() {
    HG38_SUB_RATE=hg38/simulation/substitution_rate
    SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )
    PIS=( 99 98 97 96 95 94 93 92 91 90 )

    N=12
    for i in ${!SUB_RATES[@]};
    do
        filename=${SUB_RATES[i]}
        for j in {1..100};
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_SUB_RATE}/experiment_$j/${filename}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_SUB_RATE}/experiment_$j/${filename} \
                -i ${PIS[i]} \
        &
        done
    done
}

pan_trog_benchmark && wait;
gorilla_benchmark && wait;
error_benchmark && wait;
query_length_benchmark;

