function make_pacbio_reads() {
    local pacbio=https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/PacBio_SequelII_CCS_11kb/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam
    for i in {1..10}; 
    do
        python3 scripts/make_benchmark_seqs.py \
            -a scripts/hg38_acn_conversion.txt \
            -i scripts/hg38_refseq_acn.txt \
            -o seqs/hg38/reads/pacbio/experiment_${i}_sample_reads.fasta \
            -b $pacbio \
            -x ~/Research/seqs/hg38.fna.fai \
            -n 4;
    done
};

function make_nanopore_reads() {
    local nanopore=ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam
    for i in {7..7}; 
    do
        python3 scripts/make_benchmark_seqs.py \
            -a scripts/hg38_acn_conversion.txt \
            -i scripts/hg38_refseq_acn.txt \
            -o seqs/hg38/reads/nanopore/experiment_${i}_sample_reads.fasta \
            -b $nanopore \
            -x ~/Research/seqs/hg38.fna.fai \
            -n 4;
    done
}


#make_pacbio_reads;
make_nanopore_reads
