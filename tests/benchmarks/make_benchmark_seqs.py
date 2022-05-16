from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import random
import subprocess
import pysam
import json
import pandas as pd
import argparse
import logging


def get_pysam_record(fasta_path, identifier, start, end):
    ref = pysam.FastaFile(fasta_path)
    seq = ref.fetch(identifier, start, end)
    record_id = '{0}:{1}-{2}'.format(identifier, start, end)
    fasta_record = SeqRecord(Seq(seq.upper()), record_id)
    return fasta_record



def get_fasta_record(fasta_handle, identifier, start, end):
    for record in SeqIO.parse(fasta_handle, 'fasta'):
        if identifier == record.id:
            sequence = str(record.seq[start:end])
            fasta_record = SeqRecord(Seq(sequence), record.id)
            return fasta_record


def get_ncbi_record(identifier, start, end):
    '''Retrieves seq record from NCBI'''
    Entrez.email = 'shihabdider@berkeley.edu'
    fasta_handle = Entrez.efetch(db='nucleotide', id=identifier, 
                                 rettype='fasta', retmode='text', 
                                 seq_start=start, seq_stop=end)

    fasta_record = SeqIO.read(fasta_handle, "fasta")
    fasta_handle.close()

    return fasta_record


def make_random_fasta_interval(faidx, identifier, interval_length):
    ref_lengths = {}
    with open(faidx, 'r') as handle:
        for line in handle:
            split_line = line.split('\t')
            name = split_line[0]
            length = int(split_line[1])
            ref_lengths[name] = length

    ref_length = ref_lengths[identifier]
    # get the interval
    interval_start = random.randint(0, ref_length - interval_length)
    interval_end = interval_start + interval_length - 1

    return interval_start, interval_end


def make_random_interval(identifier, interval_length):
    '''Generates random start-end for retriving from NCBI'''

    # get the ref length
    Entrez.email = 'shihabdider@berkeley.edu'
    gb_handle = Entrez.efetch(db='nucleotide', id=identifier,
                              rettype='gb', retmode='text')
    gb_record = SeqIO.read(gb_handle, "genbank")
    gb_handle.close()
    ref_length = len(gb_record.seq)

    # get the interval
    interval_start = random.randint(0, ref_length - interval_length)
    interval_end = interval_start + interval_length - 1

    return interval_start, interval_end


def records_to_fasta(records, output):
    with open(output, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        print('{} seq records saved to {}'.format(len(records), output))


def main():
    parser = argparse.ArgumentParser(
        description="Create sequence FASTA for running benchmark.s"
    )
    parser.add_argument(
        "-f", "--fasta", type=str, 
        help="path to fasta file", 
        required=True
    )
    parser.add_argument(
        "-i", "--identifiers", type=str, 
        help="path to file containing NCBI identifiers (one per line)", 
        required=True
    )
    parser.add_argument(
        "-l", "--length", type=float, 
        help=("length of the sequences"),
        default=10000
    )
    parser.add_argument(
        "-n", "--num", type=int, 
        help="number of sequences to generate per identifier", 
        default=1
    )
    parser.add_argument(
        "-o", "--output", type=str, 
        help="output path to save FASTA file", 
        required=True
    )
    args = parser.parse_args()
    faidx = args.fasta + '.fai'

    identifiers = []
    with open(args.identifiers, 'r') as handle:
        for line in handle:
            identifiers.append(line.rstrip())

    print('Retrieving random records...')
    output_records = []
    with open(args.fasta, 'r') as handle:
        for identifier in identifiers:
            for i in range(args.num):
                rand_start, rand_end = make_random_fasta_interval(
                    faidx, identifier, args.length)
                random_record = get_pysam_record(args.fasta, identifier,
                                                 rand_start, rand_end)
                print(identifier, i, random_record)
                output_records.append(random_record)

    records_to_fasta(output_records, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
