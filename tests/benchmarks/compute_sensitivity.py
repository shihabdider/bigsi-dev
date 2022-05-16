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

def compute_sensitivity(mappings, reads_chr_format):
    '''Given a list of mappings computes the sensitivity of bigsi search'''
    true_positives = 0
    total = 0
    acn_convert_df = pd.read_csv('./hg38_acn_conversion.txt', sep='\t', 
                                 header=0)
    num_buckets_hit = []
    for mapping in mappings:
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

def main():
    parser = argparse.ArgumentParser(
        description="Compute metrics on benchmark results.")
    parser.add_argument(
        "-q", "--query", type=str, 
        help="FASTA file with benchmark sequences", 
        required=True
    )
    parser.add_argument(
        "-c", "--config", type=str, 
        help="JSON file for specifying BIGSI config", 
        required=True
    )
    parser.add_argument(
        "-o", "--output", type=str, 
        help="Base file path for outputting benchmark files", 
        required=True
    )
    args = parser.parse_args()

    config = json.load(args.config)
    query_records = [record for record in SeqIO.parse(args.query, 'fasta')]

    bigsi_results = {}
    for record in query_records:
        bigsi_output = run_bigsi_query(record.seq, config)
        bigsi_results[record.id] = bigsi_output

    bigsi_results_path = args.output + '.bigsi.json'
    write_to_json(bigsi_results, bigsi_results_path)

    mashmap_results_path = args.output + '.mashmap.out'
    run_mashmap(args.query, config, output=mashmap_results_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
