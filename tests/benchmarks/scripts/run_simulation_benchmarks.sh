# Global parameters
QUERY_LENGTHS=( 1000 2000 3000 4000 5000 7500 10000 12500 15000 17500 20000 ) # 40000 80000 160000 200000 250000 300000 )
SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )

function mammal_benchmark() {
    echo "Running $1 benchmark";
    mammal_dir=$1/simulation/query_length

    for j in $( seq 1 $num_experiments );
    do
        mkdir -p outputs/${mammal_dir}/experiment_$j/;
    done

    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in $( seq 1 $num_experiments );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${mammal_dir}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${mammal_dir}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 3 \
                -m $mashmap_flag \
        &
        done
    done
}

function query_length_benchmark() {
    echo "Running query length benchmark";
    HG38_QUERY_LEN=hg38/simulation/query_length


    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in $( seq 1 $num_experiments );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 0 \
                -m $mashmap_flag \
        &
        done
    done
}

function error_and_length_benchmark() {
    echo "Running variable error rate x query length benchmark";
    HG38_SUB_RATE=hg38/simulation/error_and_query_length

    local sub_rates=( 001 002 003 004 005 )

    for j in $( seq 1 $num_experiments );
    do
        mkdir -p outputs/${HG38_SUB_RATE}/experiment_${j};
    done

    N=12
    for i in ${!QUERY_LENGTHS[@]}; do
        for j in ${!sub_rates[@]};
        do
            filename=${QUERY_LENGTHS[i]}_${sub_rates[j]}
            for k in $( seq 1 $num_experiments );
            do
                ((b=b%N)); ((b++==0)) && wait
                time python3 scripts/benchmark.py \
                    -q seqs/${HG38_SUB_RATE}/experiment_${k}/${filename}.fasta \
                    -c scripts/hg38.office.config.json \
                    -o outputs/${HG38_SUB_RATE}/experiment_${k}/${filename} \
                    -i ${sub_rates[j]} \
                    -m $mashmap_flag \
            &
            done
        done
    done
}
function error_benchmark() {
    echo "Running error rate benchmark";
    HG38_SUB_RATE=hg38/simulation/substitution_rate

    N=12
    for i in ${!SUB_RATES[@]};
    do
        filename=${SUB_RATES[i]}
        for j in $( seq 1 $num_experiments );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_SUB_RATE}/experiment_$j/${filename}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_SUB_RATE}/experiment_$j/${filename} \
                -i ${SUB_RATES[i]} \
                -m $mashmap_flag \
        &
        done
    done
}


