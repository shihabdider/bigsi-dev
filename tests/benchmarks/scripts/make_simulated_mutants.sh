SUB_RATES=( 001 002 003 004 005 006 007 008 009 010 )
NUM_EXPERIMENTS=100

function make_mutants() {
    local ref=$1
    local N=12
    for i in ${!SUB_RATES[@]};
    do
        for j in $( seq 1 $NUM_EXPERIMENTS ); 
        do
                ((job=job%N)); ((job++==0)) && wait
            python scripts/mutate_seqs.py \
                -i seqs/${1}/simulation/substitution_rate/experiment_$j/${1}_sample.fasta \
                -r ${SUB_RATES[i]} \
                -o seqs/${1}/simulation/substitution_rate/experiment_$j/${SUB_RATES[i]}.fasta \
            &
        done
    done
}

make_mutants "hg38"
#make_mutants "synth_300M"

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
