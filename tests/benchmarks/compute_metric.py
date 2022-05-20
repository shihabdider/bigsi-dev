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
import numpy as np


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


def compute_specificity(bigsi_results, mashmap_results):
    '''Computes specificity = TN/(TN+FP)
    Params: bigsi_results and mashmap_results are dicts with the same key 
    format

    A true negative is when a mashmap mapping for a given query sequence does
    not fall in a bin that is not reported by BIGSI. A false positive is when a
    query's mashmap mapping does not match the bucket reported by BIGSI.
    '''
    total_num_bins = 16
    true_negatives = 0
    false_positives = 0
    for query in mashmap_results:
        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query]

        num_matches = 0
        num_no_matches = 0 
        for mashmap_mapping in mashmap_mappings:
            if is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
                num_matches += 1
            else:
                num_no_matches += 1

        false_positives += len(bigsi_mappings) - num_matches
        true_negatives += total_num_bins - len(bigsi_mappings) - num_no_matches

    specificity = true_negatives / (true_negatives + false_positives)
    return specificity


def compute_sensitivity(bigsi_results, mashmap_results):
    '''Computes sensitivity = TP/(TP+FN)
    Params: bigsi_results and mashmap_results are dicts with the same key 
    format

    A true positive is when a mashmap mapping for a given query sequence is
    found within one of the buckets returned by the BIGSI query. A false 
    negative is when a mashmap mapping is not found in any of the buckets 
    returned by the BIGSI query.
    '''
    true_positives = 0
    false_negatives = 0
    for query in mashmap_results:
        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query] 

        for mashmap_mapping in mashmap_mappings:
            if is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
                true_positives += 1
            else:
                false_negatives += 1

    sensitivity = true_positives / (true_positives + false_negatives)
    return sensitivity


def compute_accuracy(bigsi_results, mashmap_results):
    '''Computes accuracy = TP/Total using mashmap results as truth
    set.

    Params: bigsi_results and mashmap_results are dicts with the same key
    format

    A true positive is when a mashmap mapping for a given query sequence is
    found within one of the buckets returned by the BIGSI query.
    '''

    true_positives = 0
    total = 0

    for query in mashmap_results:
        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query]

        for mashmap_mapping in mashmap_mappings:
            total += 1
            if is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
                true_positives += 1

    accuracy = true_positives/total
    return accuracy


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
        choices=['accuracy', 'sensitivity', 'specificity'],
        help="Specifies which metric to compute", 
        required=True
    )
    args = parser.parse_args()

    bigsi_results = {}
    with open(args.bigsi, 'r') as handle:
        bigsi_results = json.load(handle)
    mashmap_results = load_mashmap(args.mashmap)
    #print(len(bigsi_results.keys()), len(mashmap_results.keys()))

    if args.metric == 'accuracy':
        accuracy = compute_accuracy(bigsi_results, mashmap_results)
        print(args.bigsi, accuracy)
    elif args.metric == 'sensitivity':
        sensitivity = compute_sensitivity(bigsi_results, mashmap_results)
        print(args.bigsi, sensitivity)
    elif args.metric == 'specificity':
        specificity = compute_specificity(bigsi_results, mashmap_results)
        print(args.bigsi, specificity)
    else:
        print('Invalid metric')


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
