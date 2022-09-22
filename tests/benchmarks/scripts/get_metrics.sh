# Todo: 
# - add QUERY_LENGTHS to nanopore benchmark

function get_metric() {
    local benchmark_dir=$1
    local metric=$2
    shift 2
    local parameters=("$@")
    N=12
    for param in "${parameters[@]}"; do
        for i in $( seq 1 $num_experiments );
        do
            ((b=b%N)); ((b++==0)) && wait
            local bigsi=outputs/${benchmark_dir}/experiment_$i/${param}.bigsi.json;
            local mashmap=outputs/${benchmark_dir}/experiment_$i/${param}.mashmap.out;
            python3 scripts/compute_metric.py \
                -b ${bigsi} -m ${mashmap}  -t ${metric} -c ${bin_map} \
        &
        done
    done
}

function error_and_length_metrics() {
    local benchmark_dir=hg38/simulation/error_and_query_length
    local params=( )
    for err in ${SUB_RATES[@]}; do
        for length in ${QUERY_LENGTHS[@]}; do
            param="${length}_${err}"
            params+=( "$param" )
        done
    done

    local output=$1
    get_metric ${benchmark_dir} sensitivity ${params[@]} > metrics/${output}_sensitivity.txt
    get_metric ${benchmark_dir} specificity ${params[@]} > metrics/${output}_specificity.txt
}

function sub_rate_metrics() {
    local benchmark_dir=hg38/simulation/substitution_rate
    local output=$1
    get_metric ${benchmark_dir} sensitivity ${SUB_RATES[@]} > metrics/${output}_sensitivity.txt
    get_metric ${benchmark_dir} specificity ${SUB_RATES[@]} > metrics/${output}_specificity.txt
}

function query_length_metrics() {
    local benchmark_dir=hg38/simulation/query_length
    local output=$1
    get_metric ${benchmark_dir} sensitivity ${QUERY_LENGTHS[@]} > metrics/${output}_sensitivity.txt
    get_metric ${benchmark_dir} specificity ${QUERY_LENGTHS[@]} > metrics/${output}_specificity.txt
}

function mammal_metrics() {
    local benchmark_dir=$1/simulation/query_length
    local output=$1
    get_metric ${benchmark_dir} sensitivity ${QUERY_LENGTHS[@]} > metrics/${output}_sensitivity.txt
    get_metric ${benchmark_dir} specificity ${QUERY_LENGTHS[@]} > metrics/${output}_specificity.txt
}

function nanopore_read_metrics() {
    local METRIC=$1
    nanopore=hg38/reads/nanopore
    for i in $( seq 1 $num_experiments );
    do
            python3 scripts/compute_metric.py -b \
            outputs/${nanopore}/experiment_$i.bigsi.json -m \
            outputs/${nanopore}/experiment_$i.mashmap.out -t $METRIC \
            -c ${bin_map}
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
            -c ${bin_map}
    done
}

