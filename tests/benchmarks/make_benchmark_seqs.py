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


def get_ncbi_record(identifier, start, end):
    '''Retrieves seq record from NCBI'''
    Entrez.email = 'shihabdider@berkeley.edu'
    fasta_handle = Entrez.efetch(db='nucleotide', id=identifier, 
                                 rettype='fasta', retmode='text', 
                                 seq_start=start, seq_stop=end)

    fasta_record = SeqIO.read(fasta_handle, "fasta")
    fasta_handle.close()

    return fasta_record


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

    identifiers = []
    with open(args.identifiers, 'r') as handle:
        for line in handle:
            identifiers.append(line.rstrip())

    print('Retrieving random records...')
    output_records = []
    for identifier in identifiers:
        for i in range(args.num):
            rand_start, rand_end = make_random_interval(identifier,
                                                        args.length)
            random_record = get_ncbi_record(identifier, rand_start, rand_end)
            print(identifier, i, random_record)
            output_records.append(random_record)

    records_to_fasta(output_records, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
