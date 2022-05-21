HG38_QUERY_LEN=hg38/simulation/query_length
QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )
MIN_LENGTHS=( 500 1000 1500 2000 2500 5000 10000 20000 40000 80000 100000 125000 150000 )

N=12
for i in ${!QUERY_LENGTHS[@]};
do
    for j in {1..100};
    do
        ((b=b%N)); ((b++==0)) && wait
        time python scripts/benchmark.py \
            -q seqs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]}.fasta \
            -c scripts/hg38.office.config.json \
            -o outputs/${HG38_QUERY_LEN}/experiment_$j/${QUERY_LENGTHS[i]} \
            -i 100 \
	    -l ${MIN_LENGTHS[i]} \
	&
    done
done

#for i in {1..20}; do
#	time python benchmark.py -q seqs/hg38.experiment.$i.1000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.1000 -i 100 -l 500
#	time python benchmark.py -q seqs/hg38.experiment.$i.2000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.2000 -i 100 -l 1000
#	time python benchmark.py -q seqs/hg38.experiment.$i.3000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.3000 -i 100 -l 1500
#	time python benchmark.py -q seqs/hg38.experiment.$i.4000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.4000 -i 100 -l 2000
#	#time python benchmark.py -q seqs/hg38.experiment.$i.5000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.5000 -i 100 -l 2500
#	#time python benchmark.py -q seqs/hg38.experiment.$i.7000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.7000 -i 100 -l 3500
#	#time python benchmark.py -q seqs/hg38.experiment.$i.10000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.10000 -i 100 -l 5000
#	#time python benchmark.py -q seqs/hg38.experiment.$i.20000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.20000 -i 100 -l 10000
#	#time python benchmark.py -q seqs/hg38.experiment.$i.40000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.40000 -i 100 -l 20000
#	#time python benchmark.py -q seqs/hg38.experiment.$i.80000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.80000 -i 100 -l 40000
#	#time python benchmark.py -q seqs/hg38.experiment.$i.160000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.160000 -i 100 -l 80000 
#	#time python benchmark.py -q seqs/hg38.experiment.$i.200000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.200000 -i 100 -l 10000 
#	#time python benchmark.py -q seqs/hg38.experiment.$i.250000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.250000 -i 100 -l 125000 
#	#time python benchmark.py -q seqs/hg38.experiment.$i.300000.fasta -c hg38.office.config.json -o output/hg38.experiment.$i.300000 -i 100 -l 150000
#done

#time python benchmark.py -q seqs/synthetic_seq_300M.1000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.1000 -i 100 -l 500
#time python benchmark.py -q seqs/synthetic_seq_300M.2000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.2000 -i 100 -l 1000
#time python benchmark.py -q seqs/synthetic_seq_300M.3000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.3000 -i 100 -l 1500
#time python benchmark.py -q seqs/synthetic_seq_300M.4000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.4000 -i 100 -l 2000
#time python benchmark.py -q seqs/synthetic_seq_300M.5000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.5000 -i 100 -l 2500
#time python benchmark.py -q seqs/synthetic_seq_300M.10000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.10000 -i 100 -l 5000
#time python benchmark.py -q seqs/synthetic_seq_300M.20000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.20000 -i 100 -l 10000
#time python benchmark.py -q seqs/synthetic_seq_300M.40000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.40000 -i 100 -l 20000
#time python benchmark.py -q seqs/synthetic_seq_300M.80000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.80000 -i 100 -l 40000
#time python benchmark.py -q seqs/synthetic_seq_300M.160000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.160000 -i 100 -l 80000
#time python benchmark.py -q seqs/synthetic_seq_300M.200000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.200000 -i 100 -l 100000
#time python benchmark.py -q seqs/synthetic_seq_300M.250000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.250000 -i 100 -l 125000
#time python benchmark.py -q seqs/synthetic_seq_300M.300000.fasta -c synth_seq.office.config.json -o output/synthetic_seq_300M.300000 -i 100 -l 150000
