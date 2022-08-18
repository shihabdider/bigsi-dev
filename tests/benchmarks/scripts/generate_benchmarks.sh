#!/bin/bash

num_experiments=10;
query_config="/Users/shihabdider/Research/bigsi-dev/bigsis/hg38_32M_bins_query_config.json";
mashmap_flag=0;
script_dir="/Users/shihabdider/Research/bigsi-dev/tests/benchmarks/scripts";
commit_msg="add metrics for ${num_experiments} experiments";

echo "Computing metrics for $num_experiments experiments on $query_config BIGSI...";
echo "Running mashmap query: $mashmap_flag";

# Load functions from scripts
source "${script_dir}/run_simulation_benchmarks.sh";
source "${script_dir}/run_read_benchmarks.sh";
source "${script_dir}/get_metrics.sh";

# Simulation benchmarks

error_benchmark && wait;
query_length_benchmark && wait;
sub_rate_metrics "sub_rate_99999995_manual_adaptive_32M";
query_length_metrics "query_length_99999995_manual_adaptive_32M";

#mammal_benchmark "pan_trog" && wait;
#mammal_benchmark "gorilla" && wait;
#mammal_metrics "pan_trog";
#mammal_metrics "gorilla";

#pacbio_benchmark && wait;
#nanopore_benchmark;
#nanopore_read_metrics sensitivity > metrics/nanopore_read_sensitivities_w50.txt;
#nanopore_read_metrics specificity > metrics/nanopore_read_specificities_w50.txt; 
#pacbio_read_metrics sensitivity > metrics/pacbio_read_sensitivities_w50.txt;
#pacbio_read_metrics specificity > metrics/pacbio_read_specificities_w50.txt

# Push to remote repo
git add . && git commit -m "${commit_msg}";
git push;

echo -e "\a";
