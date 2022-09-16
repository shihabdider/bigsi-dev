import json
import argparse
import logging
import numpy as np

class Interval:
    def __init__(self, ref, start, end, bin_map):
        self.ref = ref # str
        self.start = start # int
        self.end = end # int
        self.containing_bins = self.get_containing_bins(bin_map)
        self.non_containing_bins = self.get_non_containing_bins(bin_map)

    def get_containing_bins(self, bin_map):
        '''Returns a list of the bins which contain this interval'''
        containing_bins = []
        for bin_num in bin_map:
            bin_ref = bin_map[bin_num]["refName"]
            bin_start = bin_map[bin_num]["bucketStart"]
            bin_end = bin_map[bin_num]["bucketEnd"]
            if (self.ref == bin_ref and
                int(self.start) >= int(bin_start) and
                int(self.end) <= int(bin_end)):
                containing_bins.append(bin_num)

        return containing_bins

    def get_non_containing_bins(self, bin_map):
        '''Returns a list of the bins which do not contain this interval'''
        non_containing_bins = []
        for bin_num in bin_map:
            bin_ref = bin_map[bin_num]["refName"]
            bin_start = bin_map[bin_num]["bucketStart"]
            bin_end = bin_map[bin_num]["bucketEnd"]
            if not (self.ref == bin_ref and
                int(self.start) >= int(bin_start) and
                int(self.end) <= int(bin_end)):
                non_containing_bins.append(bin_num)

        return non_containing_bins


def is_in_bigsi_bin(mashmap_mapping, bigsi_mappings):
    '''Checks if a mashmap mapping is in any of the reported bigsi bins'''
    for bigsi_mapping in bigsi_mappings:
        for mashmap_bin in mashmap_mapping.containing_bins:
            if mashmap_bin in bigsi_mapping.containing_bins:
                return True

def compute_specificity(bigsi_results, mashmap_results):
    '''Computes specificity = TN/(TN+FP)
    Params: mashmap_mappings and uigsi_mappings are dicts with the same key 
    format

    A true negative is when BIGSI correctly does not report the bins in which
    the mashmap mapping is not found (i.e intersection of the non-containing 
    bins). A false positive is when a query's mashmap mapping does not match 
    the bucket reported by BIGSI.
    '''

    true_negatives = 0
    false_positives = 0
    for query in mashmap_results:
        mashmap_containing_bins = set()
        mashmap_non_containing_bins = set()

        bigsi_containing_bins = set()
        bigsi_non_containing_bins = set()

        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query]
        for mashmap_mapping in mashmap_mappings:
           mashmap_containing_bins.update(mashmap_mapping.containing_bins)
           mashmap_non_containing_bins.update(mashmap_mapping.non_containing_bins)

        for bigsi_mapping in bigsi_mappings:
           bigsi_containing_bins.update(bigsi_mapping.containing_bins)
           bigsi_non_containing_bins.update(bigsi_mapping.non_containing_bins)

        # Remove the containing bins from non-containing bins
        mashmap_non_containing_bins = mashmap_non_containing_bins - mashmap_containing_bins
        bigsi_non_containing_bins = bigsi_non_containing_bins - bigsi_containing_bins

        query_tn = len(mashmap_non_containing_bins.intersection(bigsi_non_containing_bins))
        true_negatives += query_tn

        query_fp = len(bigsi_containing_bins - mashmap_containing_bins)
        false_positives += query_fp

    specificity = true_negatives / (true_negatives + false_positives)
    return specificity


def compute_average_specificity(bigsi_results, mashmap_results, total_num_bins):
    '''Computes average specificity across all classes
    Params: bigsi_results and mashmap_results are dicts with the same key 
    format

    A true negative is when a mashmap mapping for a given query sequence does
    not fall in a bin that is not reported by BIGSI. A false positive is when a
    query's mashmap mapping does not match the bucket reported by BIGSI.
    '''

    specificities = []
    for query in mashmap_results:
        specificity = compute_specificity(mashmap_mappings, bigsi_mappings)
        specificities.append(specificity)

    return np.mean(np.array(specificities))

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
        mashmap_containing_bins = set()

        bigsi_containing_bins = set()
        bigsi_non_containing_bins = set()

        mashmap_mappings = mashmap_results[query]
        bigsi_mappings = bigsi_results[query]
        for mashmap_mapping in mashmap_mappings:
           mashmap_containing_bins.update(mashmap_mapping.containing_bins)

        for bigsi_mapping in bigsi_mappings:
           bigsi_containing_bins.update(bigsi_mapping.containing_bins)
           bigsi_non_containing_bins.update(bigsi_mapping.non_containing_bins)

        # Remove the containing bins from non-containing bins
        bigsi_non_containing_bins = bigsi_non_containing_bins - bigsi_containing_bins

        query_tp = len(mashmap_containing_bins.intersection(bigsi_containing_bins))
        true_positives += query_tp

        query_fn = len(bigsi_non_containing_bins.intersection(mashmap_containing_bins))
        print(bigsi_containing_bins)
        print(mashmap_containing_bins)
        print(query_fn)
        false_negatives += query_fn

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

def load_bigsi_results(bigsi_output, bin_map):
    bigsi_results = {}
    with open(bigsi_output, 'r') as handle:
        bigsi_results = json.load(handle)

    for query in bigsi_results:
        bin_intervals = []
        for interval_str in bigsi_results[query]:
            if interval_str:
                split_line = interval_str.split(' ')
                bin_ref = split_line[0]
                bin_start = split_line[1]
                bin_end = split_line[2]
                bin_interval = Interval(
                        bin_ref, int(bin_start), int(bin_end),
                        bin_map)
                bin_intervals.append(bin_interval)
        bigsi_results[query] = bin_intervals

    return bigsi_results


def load_mashmap(mashmap_output, bin_map):
    mashmap_results = {}
    with open(mashmap_output, 'r') as handle:
        for line in handle:
            split_line = line.split(' ')
            record_id = split_line[0]
            mapping_ref = split_line[-5]
            mapping_start = split_line[-3]
            mapping_end = split_line[-2]
            mapping = Interval(
                    mapping_ref, int(mapping_start), int(mapping_end),
                    bin_map)
            if record_id in mashmap_results:
                mashmap_results[record_id].append(mapping)
            else:
                mashmap_results[record_id] = [mapping]

    return mashmap_results


def main():
    parser = argparse.ArgumentParser(
        description="Compute metrics on benchmark results."
    )
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
        "-c", "--config", type=str,
        help="BIGSI configuration file containing BIGSI bin map",
        required=True
    )
    parser.add_argument(
        "-t", "--metric", type=str,
        choices=['accuracy', 'sensitivity', 'specificity'],
        help="Specifies which metric to compute",
        required=True
    )

    args = parser.parse_args()
    bin_map = {}
    with open(args.config, 'r') as handle:
        bin_map = json.load(handle)

    bigsi_results = load_bigsi_results(args.bigsi, bin_map)
    mashmap_results = load_mashmap(args.mashmap, bin_map)

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
