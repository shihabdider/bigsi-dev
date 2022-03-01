'''Generates sequences from hg38 chr1 with random substitutions'''

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import random

HG38_CHR1_PATH = '/Users/shihabdider/Research/seqs/human/Homo_sapiens.GRCh38.dna.chromosome.1.fa'
#hg38_chr1 = str([record.seq for record in SeqIO.parse(HG38_CHR1_PATH, 'fasta')
#                if len(record.seq) > 100e6][0])

CHIMP_CHR1_PATH = '/Users/shihabdider/Research/seqs/mammals/chr1.fa'
#chimp_chr1 = str([record.seq for record in SeqIO.parse(CHIMP_CHR1_PATH, 'fasta')
#                if len(record.seq) > 100e6][0])

GORILLA_CHR1_PATH = '/Users/shihabdider/Research/seqs/mammals/gorilla_chr1.fasta'
#gorilla_chr1 = str([record.seq for record in SeqIO.parse(GORILLA_CHR1_PATH, 'fasta')
#                if len(record.seq) > 100e6][0])

COW_CHR1_PATH = '/Users/shihabdider/Research/seqs/mammals/cow_chr1.fasta'
#cow_chr1 = str([record.seq for record in SeqIO.parse(COW_CHR1_PATH, 'fasta')
#                if len(record.seq) > 100e6][0])

PIG_CHR1_PATH = '/Users/shihabdider/Research/seqs/mammals/pig_chr1.fasta'
#pig_chr1 = str([record.seq for record in SeqIO.parse(PIG_CHR1_PATH, 'fasta')
#                if len(record.seq) > 100e6][0])

MOUSE_CHR1_PATH = '/Users/shihabdider/Research/seqs/mammals/mouse_chr1.fasta'
mouse_chr1 = str([record.seq for record in SeqIO.parse(MOUSE_CHR1_PATH, 'fasta')
                if len(record.seq) > 100e6][0])

def generate_random_query(ref, ref_name, query_len, error_rate):
    query_start = random.randint(0, len(ref) - query_len)
    query_end = query_start + query_len

    query_seq = ref[query_start:query_end]
    random_query = []
    num_subs = 0
    for pos in range(query_len):
        if random.random() < error_rate:
            num_subs += 1
            random_nuc = random.choice(['A', 'C', 'G', 'T'])
            random_query.append(random_nuc)
        else:
            random_query.append(query_seq[pos])

    true_error_rate = num_subs/query_len
    random_query_seq = ''.join(random_query)
    random_query_name = '{}:{}_{}:{}'.format(ref_name, query_start, query_end,
                                             true_error_rate)

    random_query_rec = {'seq': random_query_seq, 'name': random_query_name}
    return random_query_rec


def main():
    num_generated_seqs = 20
    query_len = int(300e3)
    error_rate = 0.0
    for i in range(num_generated_seqs):
        random_query_rec = generate_random_query(mouse_chr1, 'mouse_chr1',
                                                 query_len,
                                                 error_rate)
        filename = str(i) + random_query_rec['name']
        random_seq_record = SeqRecord(Seq(random_query_rec['seq']),
                                      id=filename)
        path = 'tests/random_seqs/{}/{}.fasta'.format(error_rate, filename)
        with open(path, "w") as output_handle:
            SeqIO.write(random_seq_record, output_handle, "fasta")
            print(filename + ' written to ' + path)


main()
