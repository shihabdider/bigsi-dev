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
import math


def mutate_sequence(sequence, mutation_rate):
    '''Adds random base substitutions to a sequence'''

    seq_length = len(sequence)
    num_mutations = math.floor(seq_length * mutation_rate)
    mutant_pos = random.sample(range(seq_length), num_mutations)

    mutated_seq_arr = list(sequence)
    for pos in mutant_pos:
        random_nucs = ['A', 'C', 'G', 'T']
        if sequence[pos] in random_nucs:
            random_nucs.remove(sequence[pos])
            mutant_nuc = random.choice(random_nucs)
            mutated_seq_arr[pos] = mutant_nuc

    mutated_seq = ''.join(mutated_seq_arr)

    return mutated_seq


def mutate_records(records, mutation_rate):
    mutated_records = []
    for record in records:
        mutated_seq = mutate_sequence(record.seq, mutation_rate)
        mutated_record = SeqRecord(Seq(mutated_seq), id=record.id)
        mutated_records.append(mutated_record)

    return mutated_records


def records_to_fasta(records, output):
    with open(output, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        print('{} seq records saved to {}'.format(len(records), output))


def main():
    parser = argparse.ArgumentParser(
        description="Create sequence FASTA for running benchmarks"
    )
    parser.add_argument(
        "-i", "--input", type=str, 
        help=("path to FASTA file of sequences to mutate"),
        required=True
    )
    parser.add_argument(
        "-r", "--rate", type=float, 
        help=("rate to mutate sequence (from 0.01 to 0.05)"),
        required=True
    )
    parser.add_argument(
        "-o", "--output", type=str, 
        help="output path to save FASTA file", 
        required=True
    )
    args = parser.parse_args()

    output_records = []
    if args.rate <= 0:
        logging.error(
            "Mutation rate must be greater than 0"
        )
        exit(1)
        print('Mutating records...')
    else: 
        with open(args.input, 'r') as handle:
            records = [record for record in SeqIO.parse(handle, 'fasta')]
            output_records = mutate_records(records, args.rate)

    records_to_fasta(output_records, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
