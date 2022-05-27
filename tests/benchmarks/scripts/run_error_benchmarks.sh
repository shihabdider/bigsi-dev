HG38_SUB_RATE=hg38/simulation/substitution_rate
SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )
PIS=( 99 98 97 96 95 94 93 92 91 90 )

N=12
for i in ${!SUB_RATES[@]};
do
    filename=${SUB_RATES[i]}
    for j in {1..10};
    do
        ((b=b%N)); ((b++==0)) && wait
        time python3 scripts/benchmark.py \
            -q seqs/${HG38_SUB_RATE}/experiment_$j/${filename}.fasta \
            -c scripts/hg38.office.config.json \
            -o outputs/${HG38_SUB_RATE}/experiment_$j/${filename}_adaptive_error \
            -i ${PIS[i]} \
	&
    done
done

#time python benchmark.py -q seqs/hg38.random.007.fasta -c hg38.office.config.json -o output/hg38.random.007 -i 93
#time python benchmark.py -q seqs/hg38.random.008.fasta -c hg38.office.config.json -o output/hg38.random.008 -i 92
#time python benchmark.py -q seqs/hg38.random.009.fasta -c hg38.office.config.json -o output/hg38.random.009 -i 91

#time python benchmark.py -q seqs/synthetic_seq_300M.random.001.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.001 -i 99
#time python benchmark.py -q seqs/synthetic_seq_300M.random.002.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.002 -i 98
#time python benchmark.py -q seqs/synthetic_seq_300M.random.003.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.003 -i 97
#time python benchmark.py -q seqs/synthetic_seq_300M.random.004.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.004 -i 96
#time python benchmark.py -q seqs/synthetic_seq_300M.random.005.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.005 -i 95
#time python benchmark.py -q seqs/synthetic_seq_300M.random.006.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.006 -i 94
#time python benchmark.py -q seqs/synthetic_seq_300M.random.007.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.007 -i 93
#time python benchmark.py -q seqs/synthetic_seq_300M.random.008.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.008 -i 92
#time python benchmark.py -q seqs/synthetic_seq_300M.random.009.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.009 -i 91
#time python benchmark.py -q seqs/synthetic_seq_300M.random.010.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.random.010 -i 90
