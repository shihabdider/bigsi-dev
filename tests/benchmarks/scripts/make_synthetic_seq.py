import random
import sys
from Bio import SeqIO
from Bio import Seq
from Bio import SeqRecord


def make_synthetic_genome(seq_length):
    dna = ["A", "G", "C", "T"]
    genome_records = []
    if seq_length > 100e6:
        for i in range(0, seq_length, int(100e6)):
            synthetic_id = 'synthetic_seq {}'.format(i)
            synthetic_seq = make_synthetic_seq(dna, 100e6)
            synth_record = SeqRecord(Seq(synthetic_seq), synthetic_id))
            genome_records.append(synth_record)
    else:
            synthetic_id = 'synthetic_seq'
            synthetic_seq = make_synthetic_seq(dna, seq_length)
            synth_record = SeqRecord(Seq(synthetic_seq), synthetic_id))
            genome_records.append(synth_record)

    return genome_records


def make_synthetic_seq(alphabet, seq_length):
    synthetic_seq = ''
    for i in range(1, int(seq_length+1)):
        synthetic_seq += random.choice(alphabet)
        if i % 80 == 0:
            synthetic_seq += '\n'

    return synthetic_seq


def main():
    seq_length = int(sys.argv[1])*10e6
    fasta = sys.argv[2]
    print('Generating sequence of length: {0}'.format(seq_length))
    genome = make_synthetic_genome(seq_length)
    with open(fasta, 'w') as handle:
        SeqIO.write(genome, handle, 'fasta')

main()
