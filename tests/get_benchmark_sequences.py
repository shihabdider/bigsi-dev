import Bio
from Bio import Entrez
Entrez.email = 'shihabdider@gmail.com'
record = Entrez.efetch(db = 'nucleotide', id = 'NC_051849.1', rettype = 'fasta', retmode = 'text', seq_start = 33845728, seq_stop = 33848021)
print(len(record.read()))
