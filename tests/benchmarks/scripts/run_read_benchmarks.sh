MASHMAP_FLAG=$1
function nanopore_benchmark() {
    NANOPORE_READ_DIR=hg38/reads/nanopore
    for j in {1..10};
        do
            time python3 scripts/benchmark.py \
                -q seqs/${NANOPORE_READ_DIR}/experiment_${j}_sample_reads.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${NANOPORE_READ_DIR}/experiment_${j} \
                -i 5 \
                -m $MASHMAP_FLAG \
            &
    done
}

function pacbio_benchmark() {
PACBIO_READ_DIR=hg38/reads/pacbio
    for j in {1..10};
        do
            time python3 scripts/benchmark.py \
                -q seqs/${PACBIO_READ_DIR}/experiment_${j}_sample_reads.fasta \
                -c scripts/hg38.office.config.json \
                -o outputs/${PACBIO_READ_DIR}/experiment_${j} \
                -i 5 \
                -m $MASHMAP_FLAG \
            &
    done
}

