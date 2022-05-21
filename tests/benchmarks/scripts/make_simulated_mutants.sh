SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )
PIS=( 0.01 0.02 0.03 0.04 0.05 0.06 0.07 0.08 0.09 0.1 )

N=12
for i in ${!SUB_RATES[@]};
do
	for j in {1..100}; 
	do
        	((job=job%N)); ((job++==0)) && wait
		python scripts/mutate_seqs.py \
			-i seqs/hg38/simulation/substitution_rate/experiment_$j/hg38_sample.fasta \
			-r ${PIS[i]} \
			-o seqs/hg38/simulation/substitution_rate/experiment_$j/${SUB_RATES[i]}.fasta \
		&
	done
done

#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.01 -o seqs/synthetic_seq_300M.random.001.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.02 -o seqs/synthetic_seq_300M.random.002.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.03 -o seqs/synthetic_seq_300M.random.003.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.04 -o seqs/synthetic_seq_300M.random.004.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.05 -o seqs/synthetic_seq_300M.random.005.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.06 -o seqs/synthetic_seq_300M.random.006.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.07 -o seqs/synthetic_seq_300M.random.007.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.08 -o seqs/synthetic_seq_300M.random.008.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.09 -o seqs/synthetic_seq_300M.random.009.fasta
#python mutate_seqs.py -i seqs/synthetic_seq_300M.random.fasta -r 0.10 -o seqs/synthetic_seq_300M.random.010.fasta

#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.02 -o seqs/hg38.random.002.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.03 -o seqs/hg38.random.003.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.04 -o seqs/hg38.random.004.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.05 -o seqs/hg38.random.005.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.06 -o seqs/hg38.random.006.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.07 -o seqs/hg38.random.007.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.08 -o seqs/hg38.random.008.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.09 -o seqs/hg38.random.009.fasta
#python mutate_seqs.py -i seqs/hg38.random.fasta -r 0.10 -o seqs/hg38.random.010.fasta
