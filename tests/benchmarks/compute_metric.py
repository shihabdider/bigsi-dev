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


def compute_sensitivity_species(bigsi_mappings, mashmap_mappings):
    '''Given a list of mappings and dict of mashmap mappings computes the
    sensitivity of bigsi search'''
    true_positives = 0
    total = 0

    num_buckets_hit = []
    for mapping in bigsi_mappings:
        is_true_positive = False
        mapping_list = mapping.split('\t')
        seq_key = '{0}:{1}-{2}'.format(mapping_list[1], mapping_list[2],
                                       mapping_list[3])
        if seq_key in mashmap_mappings.keys():
            mashmaps = mashmap_mappings[seq_key]
            bigsi_output = mapping_list[-1]
            if bigsi_output and mashmaps:
                for mapping in mashmaps:
                    seq_ref, seq_start, seq_end = mapping.split(' ')
                    bigsi_mappings = bigsi_output.split(',')[0:-1]
                    num_buckets_hit.append(len(bigsi_mappings))
                    for bigsi_map in bigsi_mappings:
                        bigsi_map_list = bigsi_map.split(' ')
                        bin_ref = bigsi_map_list[0]
                        bin_start = int(bigsi_map_list[1])
                        bin_end = int(bigsi_map_list[2])
                        is_positive = (
                            seq_ref == bin_ref and
                            int(seq_start) >= bin_start and
                            int(seq_end) <= bin_end
                        )

                        if is_positive:
                            is_true_positive = True

                if is_true_positive:
                    true_positives += 1

                total += 1

    sensitivity = true_positives/total
    num_multibin_hits = [num_bins for num_bins in num_buckets_hit if num_bins > 20]
    print(num_buckets_hit, len(num_multibin_hits), len(num_buckets_hit))
    return sensitivity


def compute_sensitivity(bigsi_results, mashmap_results):
    '''Computes sensitivity = TP/(TP + FN) using mashmap results as truth
    set. 

    Params: bigsi_results and mashmap_results are dicts with the same key 
    format

    A true positive is when all the mashmap mappings for a given query sequence 
    are found within the buckets returned by the BIGSI query'''

    true_positives = 0
    total = 0

    num_buckets_hit = []
    for mapping in mashmap_results:
        is_true_positive = False
        mapping_list = mapping.split('\t')
        read_ref = mapping_list[1]
        read_start = int(mapping_list[2])
        read_end = int(mapping_list[3])
        bigsi_output = mapping_list[-1]
        if bigsi_output:
            bigsi_mappings = bigsi_output.split(',')[0:-1]
            num_buckets_hit.append(len(bigsi_mappings))
            for bigsi_map in bigsi_mappings:
                bigsi_map_list = bigsi_map.split(' ')
                mapped_ref_name = acn_convert_df[
                    acn_convert_df['RefSeq-Accn'].str.contains(
                        bigsi_map_list[0])][reads_chr_format].values[0]
                bin_start = int(bigsi_map_list[1])
                bin_end = int(bigsi_map_list[2])
                is_positive = (read_ref == mapped_ref_name and
                                read_start >= bin_start and
                                read_end <= bin_end)

                if is_positive:
                    is_true_positive = True

        if is_true_positive:
            true_positives += 1

        total += 1

    sensitivity = true_positives/total
    num_multibin_hits = [num_bins for num_bins in num_buckets_hit if num_bins > 20]
    print(num_buckets_hit, len(num_multibin_hits), len(num_buckets_hit))
    return sensitivity


def load_mashmap(mashmap_output):
    mashmap_results = {}
    with open(mashmap_output, 'r') as handle:
        for line in handle:
            split_line = line.split(' ')
            record_id = split_line[0]
            mapping = split_line[-5:-1]
            if mashmap_results[record_id]:
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

    bigsi_results = json.load(args.bigsi)
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
