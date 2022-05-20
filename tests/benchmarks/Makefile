#all: figures/%.png

.PHONY: all clean

HG38_SUB_RATE_SIM=hg38/simulation/substitution_rate
HG38_QUERY_LEN_SIM=hg38/simulation/query_length

METRICS_HG38_SUB_RATE_SIM := metrics/${HG38_SUB_RATE_SIM}/trials_sensitivity.txt \
			metrics/${HG38_SUB_RATE_SIM}/trials_specificity.txt \

METRICS_HG38_QUERY_LEN_SIM := metrics/${HG38_QUERY_LEN_SIM}/trials_sensitivity.txt \
			metrics/${HG38_QUERY_LEN_SIM}/trials_specificity.txt

METRICS_HG38_SIM := ${METRICS_HG38_SUB_RATE_SIM} ${METRICS_HG38_QUERY_LEN_SIM}

figures/flashmap_accuracy_hg38_simulation_trials.png: ${METRICS_HG38_SIM}
	scripts/make_benchmark_figures.py $^ $@

${METRICS_HG38_SUB_RATE_SIM}: outputs/${HG38_SUB_RATE_SIM}/error_rate/experiment_*/%.bigsi.json outputs/${HG38_SUB_RATE_SIM}/experiment_*/%.mashmap.out
	scripts/get_metrics.sh sensitivity outputs/${HG38_SUB_RATE_SIM} > metrics/${HG38_SUB_RATE_SIM}/trials_sensitivity.txt
	scripts/get_metrics.sh specificity outputs/${HG38_SUB_RATE_SIM} > metrics/${HG38_SUB_RATE_SIM}/trials_specificity.txt

${METRICS_HG38_QUERY_LEN_SIM}: outputs/${HG38_QUERY_LEN_SIM}/error_rate/experiment_*/%.bigsi.json outputs/${HG38_QUERY_LEN_SIM}/experiment_*/%.mashmap.out
	scripts/get_metrics.sh sensitivity outputs/${HG38_QUERY_LEN_SIM} > metrics/${HG38_QUERY_LEN_SIM}/trials_sensitivity.txt
	scripts/get_metrics.sh specificity outputs/${HG38_QUERY_LEN_SIM} > metrics/${HG38_QUERY_LEN_SIM}/trials_specificity.txt
