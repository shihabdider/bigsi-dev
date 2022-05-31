QUERY_LENGTHS=( 1000 2000 3000 4000 5000 10000 20000 40000 80000 160000 200000 250000 300000 )

function pan_trog_dataset() {
    local pan_trog_dir=pan_trog/simulation
    for i in ${QUERY_LENGTHS[@]};
    do
        for j in {1..100};
        do
            mkdir -p seqs/${pan_trog_dir}/experiment_${j};
            python3 scripts/make_benchmark_seqs.py \
                -i scripts/pan_trog_acn.txt \
                -l ${i} \
                -n 4 \
                -o seqs/${pan_trog_dir}/experiment_${j}/${i}.fasta \
                -f ~/Research/seqs/pan_trog.fasta \
                -x ~/Research/seqs/pan_trog.fasta.fai
        done
    done
}


function gorilla_dataset() {
    local gorilla_dir=gorilla/simulation
    for i in ${QUERY_LENGTHS[@]};
    do
        for j in {1..100};
        do
            mkdir -p seqs/${gorilla_dir}/experiment_${j};
            python3 scripts/make_benchmark_seqs.py \
                -i scripts/gorilla_acn.txt \
                -l ${i} \
                -n 4 \
                -o seqs/${gorilla_dir}/experiment_${j}/${i}.fasta \
                -f ~/Research/seqs/gorilla.fasta \
                -x ~/Research/seqs/gorilla.fasta.fai
        done
    done
}

pan_trog_dataset;
gorilla_dataset
