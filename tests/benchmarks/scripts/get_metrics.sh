# Todo: 
# - add lengths to nanopore benchmark

function run_benchmark() {
    local benchmark_dir=$1
    local metric=$2
    shift 2
    local parameters=("$@")
    for param in "${parameters[@]}"; do
        for i in $( seq 1 $num_experiments );
        do
            local bigsi=outputs/${benchmark_dir}/experiment_$i/${param}.bigsi.json;
            local mashmap=outputs/${benchmark_dir}/experiment_$i/${param}.mashmap.out;
            python3 scripts/compute_metric.py \
                -b ${bigsi} -m ${mashmap}  -t ${metric} -c ${query_config}
        done
    done
}

function sub_rate_metrics() {
    local benchmark_dir=hg38/simulation/substitution_rate
    local errs=( 001 002 003 004 005 006 007 008 009 010 )
    local output=$1
    run_benchmark ${benchmark_dir} sensitivity ${errs[@]} > metrics/${output}_sensitivity.txt
    run_benchmark ${benchmark_dir} specificity ${errs[@]} > metrics/${output}_specificity.txt
}

function query_length_metrics() {
    local benchmark_dir=hg38/simulation/query_length
    local lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    local output=$1
    run_benchmark ${benchmark_dir} sensitivity ${lengths[@]} > metrics/${output}_sensitivity.txt
    run_benchmark ${benchmark_dir} specificity ${lengths[@]} > metrics/${output}_specificity.txt
}

function mammal_metrics() {
    local benchmark_dir=$1/simulation/
    local lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    local output=$1
    run_benchmark ${benchmark_dir} sensitivity ${lengths[@]} > metrics/${output}_sensitivity.txt
    run_benchmark ${benchmark_dir} specificity ${lengths[@]} > metrics/${output}_specificity.txt
}

function nanopore_read_metrics() {
    local METRIC=$1
    nanopore=hg38/reads/nanopore
    for i in $( seq 1 $num_experiments );
    do
            python3 scripts/compute_metric.py -b \
            outputs/${nanopore}/experiment_$i.bigsi.json -m \
            outputs/${nanopore}/experiment_$i.mashmap.out -t $METRIC \
            -c ${query_config}
    done
}

function pacbio_read_metrics() {
    local METRIC=$1
    pacbio=hg38/reads/pacbio
    for i in $( seq 1 $num_experiments );
    do
            python3 scripts/compute_metric.py -b \
            outputs/${pacbio}/experiment_$i.bigsi.json -m \
            outputs/${pacbio}/experiment_$i.mashmap.out -t $METRIC \
            -c ${query_config}
    done
}

