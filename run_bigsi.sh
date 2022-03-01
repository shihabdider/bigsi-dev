#REF="../../seqs/human/ch17.fasta"
#OUTPUT=bigsis/hg38_chr17
REF="../../seqs/human/hg38/ncbi-genomes-2021-11-16/hg38.fna"
#OUTPUT=hg38_whole_genome_inexact_hashes
#REF="../../seqs/human/Homo_sapiens.GRCh38.dna.chromosome.1.fa"
OUTPUT=bigsis/hg38_whole_genome_005.bin
###
#### Make the bigsi
#node bigsi.js -r $REF -o $OUTPUT
#
# Paths must be absolute or relative to query_bigsi.js
#QUERY="../../seqs/human/brca1.fasta"
#QUERY="../../seqs/human/ltr7.fasta"
#QUERY="../../seqs/human/l1td1.fasta"
#QUERY="../../seqs/human/MER57A-int.fasta"
#BIGSI="bigsis/hg38_whole_genome.bin"
#CONFIG="../bigsi.query.config.json"

#QUERY="../../seqs/human/ACADM_hg38chr1.fasta"
#BIGSI="bigsis/hg38_chr1_error_rate.bin"
#BIGSI="bigsis/hg38_chr1_error_rate_low_window.bin"
#BIGSI="bigsis/hg38_chr1_error_rate_005.bin"
#BIGSI="bigsis/hg38_chr1_error_rate_003.bin"
#BIGSI="bigsis/hg38_chr1_error_rate_005.bin"
BIGSI="bigsis/hg38_chr1_error_rate_005.bin"
CONFIG="../bigsi.random.query.config.json"
#QUERY="../../seqs/mammals/mouse_align.fasta"
QUERY="../../seqs/human/nanopore_lr.fasta"
#for QUERY in tests/random_seqs/0.0/*.fasta; do
## Query the bigsi
#    node bin/query_bigsi.js -q $QUERY -b $BIGSI -c $CONFIG
#done

# Query the bigsi
node bin/query_bigsi.js -q $QUERY -b $BIGSI -c $CONFIG
