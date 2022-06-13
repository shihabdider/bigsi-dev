num_experiments=$1
echo "Computing metrics for $num_experiments experiments..."

function sub_rate_metrics() {
    local METRIC=$1
    HG38_SUB_RATE=hg38/simulation/substitution_rate
    errs=( 001 002 003 004 005 006 007 008 009 010 )
    for err in "${errs[@]}"; do
        for i in $( seq 1 $num_experiments );
        do
            python3 scripts/compute_metric.py -b outputs/${HG38_SUB_RATE}/experiment_$i/${err}.bigsi.json -m outputs/${HG38_SUB_RATE}/experiment_$i/${err}.mashmap.out -t $METRIC
        done
    done
}


function query_length_metrics() {
    local METRIC=$1
    HG38_QUERY_LEN=hg38/simulation/query_length
    lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    for length in "${lengths[@]}"; do
        for i in $( seq 1 $num_experiments );
        do
            python3 scripts/compute_metric.py -b outputs/${HG38_QUERY_LEN}/experiment_$i/${length}.bigsi.json -m outputs/${HG38_QUERY_LEN}/experiment_$i/${length}.mashmap.out -t $METRIC
        done
    done
}


function nanopore_read_metrics() {
    local METRIC=$1
    nanopore=hg38/reads/nanopore
    for i in $( seq 1 $num_experiments );
    do
            python3 scripts/compute_metric.py -b \
            outputs/${nanopore}/experiment_$i.bigsi.json -m \
            outputs/${nanopore}/experiment_$i.mashmap.out -t $METRIC
    done
}

function pacbio_read_metrics() {
    local METRIC=$1
    pacbio=hg38/reads/pacbio
    for i in $( seq 1 $num_experiments );
    do
            python3 scripts/compute_metric.py -b \
            outputs/${pacbio}/experiment_$i.bigsi.json -m \
            outputs/${pacbio}/experiment_$i.mashmap.out -t $METRIC
    done
}

function pan_trog_metrics() {
    local METRIC=$1
    pan_trog_dir=pan_trog/simulation/
    lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    for length in "${lengths[@]}"; 
    do
        for i in $( seq 1 $num_experiments );
        do
            python3 scripts/compute_metric.py -b \
            outputs/${pan_trog_dir}/experiment_$i/${length}.bigsi.json -m \
            outputs/${pan_trog_dir}/experiment_$i/${length}.mashmap.out \
            -t $METRIC
        done
    done
}

function gorilla_metrics() {
    local METRIC=$1
    gorilla_dir=gorilla/simulation/
    lengths=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
    for length in "${lengths[@]}"; do
        for i in $( seq 1 $num_experiments );
        do
            python3 scripts/compute_metric.py -b \
            outputs/${gorilla_dir}/experiment_$i/${length}.bigsi.json -m \
            outputs/${gorilla_dir}/experiment_$i/${length}.mashmap.out \
            -t $METRIC
        done
    done
}

#sub_rate_metrics sensitivity > metrics/adaptive_error_error_sensitivity_90CI.txt
#sub_rate_metrics specificity > metrics/adaptive_error_error_specificity_90CI.txt
#query_length_metrics sensitivity > metrics/adaptive_error_length_sensitivity_90CI.txt
#query_length_metrics specificity > metrics/adaptive_error_length_specificity_90CI.txt

#nanopore_read_metrics sensitivity > metrics/nanopore_read_sensitivities_003.txt;
#nanopore_read_metrics specificity > metrics/nanopore_read_specificities_003.txt; 
#pacbio_read_metrics sensitivity > metrics/pacbio_read_sensitivities_003.txt;
#pacbio_read_metrics specificity > metrics/pacbio_read_specificities_003.txt

#pan_trog_metrics sensitivity > metrics/pan_trog_sensitivities_90CI.txt;
#pan_trog_metrics specificity > metrics/pan_trog_specificities_90CI.txt; 
#gorilla_metrics sensitivity > metrics/gorilla_sensitivities_03_unbound.txt;
#gorilla_metrics specificity > metrics/gorilla_specificities_03_unbound.txt

#python compute_metric.py -b output/synthetic_seq_300M.random.001.bigsi.json -m output/synthetic_seq_300M.random.001.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.002.bigsi.json -m output/synthetic_seq_300M.random.002.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.003.bigsi.json -m output/synthetic_seq_300M.random.003.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.004.bigsi.json -m output/synthetic_seq_300M.random.004.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.005.bigsi.json -m output/synthetic_seq_300M.random.005.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.006.bigsi.json -m output/synthetic_seq_300M.random.006.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.007.bigsi.json -m output/synthetic_seq_300M.random.007.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.008.bigsi.json -m output/synthetic_seq_300M.random.008.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.009.bigsi.json -m output/synthetic_seq_300M.random.009.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.random.010.bigsi.json -m output/synthetic_seq_300M.random.010.mashmap.out -t $METRIC

#python compute_metric.py -b output/synthetic_seq_300M.1000.bigsi.json -m output/synthetic_seq_300M.1000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.2000.bigsi.json -m output/synthetic_seq_300M.2000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.3000.bigsi.json -m output/synthetic_seq_300M.3000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.4000.bigsi.json -m output/synthetic_seq_300M.4000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.5000.bigsi.json -m output/synthetic_seq_300M.5000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.10000.bigsi.json -m output/synthetic_seq_300M.10000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.20000.bigsi.json -m output/synthetic_seq_300M.20000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.40000.bigsi.json -m output/synthetic_seq_300M.40000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.80000.bigsi.json -m output/synthetic_seq_300M.80000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.160000.bigsi.json -m output/synthetic_seq_300M.160000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.200000.bigsi.json -m output/synthetic_seq_300M.200000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.250000.bigsi.json -m output/synthetic_seq_300M.250000.mashmap.out -t $METRIC
#python compute_metric.py -b output/synthetic_seq_300M.300000.bigsi.json -m output/synthetic_seq_300M.300000.mashmap.out -t $METRIC
