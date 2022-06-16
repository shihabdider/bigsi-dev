NUM_EXPERIMENTS=$1
MASHMAP_FLAG=$2
echo "Running $NUM_EXPERIMENTS experiments...";
echo "Running mashmap query: $MASHMAP_FLAG"

# Global parameters
QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )

function mammal_benchmark() {
    echo "Running $1 benchmark";
    mammal_dir=$1/simulation

    for j in $( seq 1 $NUM_EXPERIMENTS );
    do
        mkdir -p outputs/${mammal_dir}/experiment_$j/;
    done

    N=12
    for i in ${!QUERY_LENGTHS[@]};
    do
        for j in $( seq 1 $NUM_EXPERIMENTS );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${mammal_dir}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${mammal_dir}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 5 \
                -m $MASHMAP_FLAG \
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
        for j in $( seq 1 $NUM_EXPERIMENTS );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]} \
                -i 0 \
                -m $MASHMAP_FLAG \
        &
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
        for j in $( seq 1 $NUM_EXPERIMENTS );
        do
            ((b=b%N)); ((b++==0)) && wait
            time python3 scripts/benchmark.py \
                -q seqs/${HG38_SUB_RATE}/experiment_$j/${filename}.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${HG38_SUB_RATE}/experiment_$j/${filename} \
                -i ${SUB_RATES[i]} \
                -m $MASHMAP_FLAG \
        &
        done
    done
}

mammal_benchmark pan_trog && wait;
mammal_benchmark gorilla && wait;
error_benchmark && wait;
query_length_benchmark;

