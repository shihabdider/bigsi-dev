#QUERY_LENGTHS=( 7500 12500 15000 17500 ) # 40000 80000 160000 200000 250000 300000 )
QUERY_LENGTHS=( 1000 2000 3000 4000 5000 7500 10000 12500 15000 17500 20000 ) # 40000 80000 160000 200000 250000 300000 )
num_experiments=100
num_queries_per_seq=4

function mammal_dataset() {
    local mammal_dir=$1/simulation/query_length;
    mkdir -p $mammal_dir;

    N=12
    for i in ${QUERY_LENGTHS[@]};
    do
        for j in $( seq 1 $num_experiments );
        do
            ((b=b%N)); ((b++==0)) && wait
            mkdir -p seqs/${mammal_dir}/experiment_${j};
            python3 scripts/make_benchmark_seqs.py \
                -i scripts/${1}_acn.txt \
                -l ${i} \
                -n $num_queries_per_seq \
                -o seqs/${mammal_dir}/experiment_${j}/${i}.fasta \
                -f ~/Research/seqs/$1.fasta \
                -x ~/Research/seqs/$1.fasta.fai \
        &
        done
    done
}

function make_reads() {
    local sequencer=$1
    local read_dir=$2
    for i in $( seq 1 $num_experiments );
    do
        python3 scripts/make_benchmark_seqs.py \
            -a scripts/hg38_acn_conversion.txt \
            -i scripts/hg38_refseq_acn.txt \
            -o seqs/hg38/reads/${sequencer}/experiment_${i}_sample_reads.fasta \
            -b $read_dir \
            -x ~/Research/seqs/hg38.fna.fai \
            -n $num_queries_per_seq;
    done
}

function make_error_rate_dataset() {
    for i in $( seq 1 $num_experiments );
    do
        python scripts/make_benchmark_seqs.py \
            -i scripts/hg38_refseq_acn.txt \
            -o seqs/hg38/simulation/substitution_rate/experiment_$i/hg38_sample.fasta \
            -f ~/Research/seqs/hg38.fna \
            -x ~/Research/seqs/hg38.fna.fai \
            -n $num_queries_per_seq;
    done
}

function make_error_and_length_dataset() {
    #local sub_rates=( 001 002 003 004 005 006 )
    local sub_rates=( 007 008 009 010 )
    local num_experiments=100

    for n in $( seq 1 $num_experiments );
    do
        mkdir -p seqs/hg38/simulation/error_and_query_length/experiment_${n};
    done

    N=12;
    for i in ${!QUERY_LENGTHS[@]}; do
        for j in ${!sub_rates[@]}; do
            local filename=${QUERY_LENGTHS[i]}_${sub_rates[j]}
            for k in $( seq 1 $num_experiments ); do
                ((b=b%N)); ((b++==0)) && wait
                python3 scripts/mutate_seqs.py \
                    -i seqs/hg38/simulation/query_length/experiment_${k}/${QUERY_LENGTHS[i]}.fasta \
                    -r ${sub_rates[j]} \
                    -o seqs/hg38/simulation/error_and_query_length/experiment_${k}/${filename}.fasta \
                &
            done
        done
    done
}

#function make_synth_dataset() {
#}

#make_error_and_length_dataset
mammal_dataset "pan_trog";
mammal_dataset "gorilla";
#mammal_dataset "dog";
#mammal_dataset "bas";
#mammal_dataset "bonobo";

#mammal_dataset "hg38";
#make_error_rate_dataset()

#pacbio_read_dir="https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/PacBio_SequelII_CCS_11kb/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam"
#nanopore_read_dir="ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"
#make_reads "pacbio" $pacbio_read_dir;
#make_reads "nanopore" $nanopore_read_dir;

#python make_benchmark_seqs.py -i synth_acn.txt -o seqs/synthetic_seq_300M.random.fasta -f seqs/synthetic_seq_300M.fasta -x seqs/synthetic_seq_300M.fasta.fai -n 1000;

