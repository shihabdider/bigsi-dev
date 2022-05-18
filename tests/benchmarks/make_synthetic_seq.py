import random

dna = ["A", "G", "C", "T"]
synthetic_seq = ''
seq_length = 3*10e8

for i in range(1, int(seq_length+1)):
    synthetic_seq += random.choice(dna)
    if i % 80 == 0:
        synthetic_seq += '\n'

#print('>synthetic_seq {0}'.format(seq_length))
#print(synthetic_seq)

fasta = 'seqs/synthetic_seq_300M.fasta'
with open(fasta, 'w') as handle:
    handle.write('>synthetic_seq {0}'.format(seq_length))
    handle.write(synthetic_seq)
