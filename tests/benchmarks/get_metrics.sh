METRIC=$1
#python compute_metric.py -b output/hg38_001.bigsi.json -m output/hg38_001.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_002.bigsi.json -m output/hg38_002.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_003.bigsi.json -m output/hg38_003.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_004.bigsi.json -m output/hg38_004.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_005.bigsi.json -m output/hg38_005.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_006.bigsi.json -m output/hg38_006.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38.random.007.bigsi.json -m output/hg38.random.007.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38.random.008.bigsi.json -m output/hg38.random.008.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38.random.009.bigsi.json -m output/hg38.random.009.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_010.bigsi.json -m output/hg38_010.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_5KB.bigsi.json -m output/hg38_5KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_7KB.bigsi.json -m output/hg38_7KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_20KB.bigsi.json -m output/hg38_20KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_40KB.bigsi.json -m output/hg38_40KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_80KB.bigsi.json -m output/hg38_80KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_160KB.bigsi.json -m output/hg38_160KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_200KB.bigsi.json -m output/hg38_200KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_250KB.bigsi.json -m output/hg38_250KB.mashmap.out -t $METRIC
#python compute_metric.py -b output/hg38_300KB.bigsi.json -m output/hg38_300KB.mashmap.out -t $METRIC
errs=( 001 002 003 004 005 006 007 008 009 010 )
for err in "${errs[@]}"; do
	echo "\n"
	for i in {1..20};
	do
		python compute_metric.py -b output/hg38.experiment.$i.$err.bigsi.json -m output/hg38.experiment.$i.$err.mashmap.out -t $METRIC
	done
done
