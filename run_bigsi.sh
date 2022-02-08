#REF="../../seqs/human/ch17.fasta"
#OUTPUT=bigsis/hg38_chr17
#REF="../../seqs/human/hg38/ncbi-genomes-2021-11-16/hg38.fna"
#OUTPUT=hg38_whole_genome_exact 

# Make the bigsi
#node bigsi.js -r $REF -o $OUTPUT; alarm

# Paths must be absolute or relative to query_bigsi.js
#QUERY="../../seqs/human/brca1.fasta"
#QUERY="../../seqs/human/ltr7.fasta"
#QUERY="../../seqs/human/l1td1.fasta"
QUERY="../../seqs/human/MER57A-int.fasta"
BIGSI="bigsis/hg38_whole_genome.bin"
CONFIG="../bigsi.query.config.json"
# Query the bigsi
node bin/query_bigsi.js -q $QUERY -b $BIGSI -c $CONFIG
