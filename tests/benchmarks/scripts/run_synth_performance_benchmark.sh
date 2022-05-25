BIGSI_DIR=~/Research/bigsi-dev/tests/benchmarks/bigsis

function benchmark_bigsi_size_time() {
	BIGSI_SIZES=( 100 400 800 1200 1600 2000 2500 3000 )
	QUERY_SIZE=300000
	METRIC_FILE=metrics/synth_bigsi_size_performance_times.txt

	rm $METRIC_FILE # delete old metric file
	for BIGSI_SIZE in ${BIGSI_SIZES[@]};
	do
		for i in {1..100};
		do
			{ gtime -f "%e" node ~/Research/bigsi-dev/bin/query_bigsi.js \
				-b $BIGSI_DIR/synth_${BIGSI_SIZE}.bigsi.bin -c $BIGSI_DIR/synth_${BIGSI_SIZE}.query.bigsi.config.json -q seqs/synthetic/performance/synth_query_${BIGSI_SIZE}_${QUERY_SIZE}.fasta; } \
				2>> $METRIC_FILE
		done
	done
};

function benchmark_query_size_time() {
	QUERY_SIZES=( 1000 2000 4000 8000 16000 32000 64000 128000 256000 300000 )
	BIGSI_SIZE=3000
	METRIC_FILE=metrics/synth_query_size_performance_times.txt

	rm $METRIC_FILE # delete old metric file
	for QUERY_SIZE in ${QUERY_SIZES[@]};
	do
		for i in {1..100};
		do
			{ gtime -f "%e" node ~/Research/bigsi-dev/bin/query_bigsi.js \
				-b ${BIGSI_DIR}/synth_${BIGSI_SIZE}.bigsi.bin -c ${BIGSI_DIR}/synth_${BIGSI_SIZE}.query.bigsi.config.json -q seqs/synthetic/performance/synth_query_${BIGSI_SIZE}_${QUERY_SIZE}.fasta; } \
				2>> $METRIC_FILE
		done
	done
};

benchmark_bigsi_size_time;
#benchmark_query_size_time
