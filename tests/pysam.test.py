import pysam

url = "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"

with pysam.AlignmentFile(url, "rb") as samfile:
    if read.query_sequence is not None and read.query_length > 5000:
        print(read.query_name, '\t', read.query_sequence, '\t', 
              read.query_length)

