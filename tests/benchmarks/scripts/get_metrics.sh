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

sub_rate_metrics sensitivity > metrics/adaptive_error_error_sensitivity_comp.txt
sub_rate_metrics specificity > metrics/adaptive_error_error_specificity_comp.txt
query_length_metrics sensitivity > metrics/adaptive_error_length_sensitivity_comp.txt
query_length_metrics specificity > metrics/adaptive_error_length_specificity_comp.txt

#nanopore_read_metrics sensitivity > metrics/nanopore_read_sensitivities_003.txt;
#nanopore_read_metrics specificity > metrics/nanopore_read_specificities_003.txt; 
#pacbio_read_metrics sensitivity > metrics/pacbio_read_sensitivities_003.txt;
#pacbio_read_metrics specificity > metrics/pacbio_read_specificities_003.txt

#pan_trog_metrics sensitivity > metrics/pan_trog_sensitivities_comp.txt;
#pan_trog_metrics specificity > metrics/pan_trog_specificities_comp.txt; 
#gorilla_metrics sensitivity > metrics/gorilla_sensitivities_comp.txt;
#gorilla_metrics specificity > metrics/gorilla_specificities_comp.txt
