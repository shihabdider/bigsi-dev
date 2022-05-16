METRIC=$1
python compute_metric.py -b output/hg38_001.bigsi.json -m output/hg38_001.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_002.bigsi.json -m output/hg38_002.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_003.bigsi.json -m output/hg38_003.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_004.bigsi.json -m output/hg38_004.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_005.bigsi.json -m output/hg38_005.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_006.bigsi.json -m output/hg38_006.mashmap.out -t $METRIC
python compute_metric.py -b output/hg38_010.bigsi.json -m output/hg38_010.mashmap.out -t $METRIC
