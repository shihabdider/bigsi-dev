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
import importlib


def is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
    '''Checks if a mashmap mapping is in any of the reported bigsi bins'''
    map_ref, map_start, map_end = mashmap_mapping.split(' ')
    for bigsi_mapping in bigsi_mappings:
        if bigsi_mapping:
            bin_ref, bin_start, bin_end = bigsi_mapping.split(' ')
            if (map_ref == bin_ref and
                int(map_start) >= int(bin_start) and
                int(map_end) <= int(bin_end)):
                return True

    return False


def compute_sensitivity(bigsi_results, mashmap_results):
    '''Computes sensitivity = TP/(TP + FN) using mashmap results as truth
    set. 

    Params: bigsi_results and mashmap_results are dicts with the same key 
    format

    A true positive is when all the mashmap mappings for a given query sequence 
    are found within the buckets returned by the BIGSI query'''

    true_positives = 0
    total = 0

    for query in mashmap_results:
        is_true_positive = True
        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query] 

        for mashmap_mapping in mashmap_mappings:
            if not is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
                is_true_positive = False

        if is_true_positive:
            true_positives += 1

        total += 1

    sensitivity = true_positives/total
    return sensitivity


def load_mashmap(mashmap_output):
    mashmap_results = {}
    with open(mashmap_output, 'r') as handle:
        for line in handle:
            split_line = line.split(' ')
            record_id = split_line[0]
            mapping_ref = split_line[-5]
            mapping_start = split_line[-3]
            mapping_end = split_line[-2]
            mapping = ' '.join([mapping_ref, mapping_start, mapping_end])
            if record_id in mashmap_results:
                mashmap_results[record_id].append(mapping)
            else:
                mashmap_results[record_id] = [mapping]

    return mashmap_results


def main():
    parser = argparse.ArgumentParser(
        description="Compute metrics on benchmark results.")
    parser.add_argument(
        "-b", "--bigsi", type=str, 
        help="JSON with BIGSI results", 
        required=True
    )
    parser.add_argument(
        "-m", "--mashmap", type=str, 
        help="Mashmap output file (truth set) corresponding to BIGSI results", 
        required=True
    )
    parser.add_argument(
        "-t", "--metric", type=str, 
        choices=['sensitivity', 'specificity'],
        help="Specifies which metric to compute", 
        required=True
    )
    args = parser.parse_args()

    bigsi_results = {}
    with open(args.bigsi, 'r') as handle:
        bigsi_results = json.load(handle)
    mashmap_results = load_mashmap(args.mashmap)

    if args.metric == 'sensitivity':
        sensitivity = compute_sensitivity(bigsi_results, mashmap_results)
        print(sensitivity)
    elif args.metric == 'specificity':
        pass
        # compute_specificity(bigsi_results, mashmap_results)
    else:
        print('Metric must be either sensitivity or specificity')


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
