'''
Complexity calculations for BIGSI and Mashmap
'''

import math
from scipy.stats import binom

def bits_to_mb(bits):
    return bits/(8*10**6)


def compute_num_minimizers(seq_length, window_size):
    num_minimizers = math.ceil(2*seq_length/window_size)
    return num_minimizers


def error_to_containment(error_rate, kmer_length):
    containment_score = math.exp(-1*error_rate*kmer_length)
    return containment_score


def containment_to_error(containment_score, kmer_length):
    return -1/kmer_length * math.log(containment_score)


def compute_num_hash_funcs(bf_size, num_inserted):
    '''Computes optimal number of hash functions for a bloom filter'''
    return math.ceil((bf_size/num_inserted)*math.log(2))


def log_binomial(n, k):
    return (math.lgamma(n+1) - math.lgamma(k+1) - math.lgamma(n-k+1))


def compute_prob_n_kmers_not_in_minimizer_set(n, seq_length, window_size):
    num_seq_kmers = seq_length - 16 + 1
    num_minimizers = compute_num_minimizers(seq_length, window_size)
    prob_n_kmers_not_in_minimizer_set = (
        log_binomial(num_seq_kmers - n, num_seq_kmers - num_minimizers - n) -
        log_binomial(num_seq_kmers, num_seq_kmers - num_minimizers)
    )

    return math.exp(prob_n_kmers_not_in_minimizer_set)


def compute_prob_n_kmers_in_minimizer_set(n, seq_length, window_size):
    num_seq_kmers = seq_length - 16 + 1
    num_minimizers = compute_num_minimizers(seq_length, window_size)
    prob_n_kmers_in_minimizer_set = (
        log_binomial(num_seq_kmers - n, num_minimizers - n) -
        log_binomial(num_seq_kmers, num_minimizers)
    )

    return math.exp(prob_n_kmers_in_minimizer_set)


def compute_winnow_false_neg(query_size, target_size,
                             window_size, error_rate):
    '''Computes the false negative probability arising because of the
    winnnowing process for a single bin'''

    num_query_minimizers = compute_num_minimizers(query_size, window_size)
    containment_thresh = error_to_containment(error_rate, kmer_length=16)
    num_containment_queries = int(math.ceil(
        containment_thresh*num_query_minimizers)
    )

    false_negative_prob = 0
    for i in range(num_query_minimizers - num_containment_queries+1, num_query_minimizers+1):
        prob_unmatched_kmers_in_query = compute_prob_n_kmers_in_minimizer_set(
            i, query_size, window_size
        )
        prob_unmatched_kmers_not_in_target = compute_prob_n_kmers_not_in_minimizer_set(
            i, target_size, window_size
        )

        false_negative_prob += prob_unmatched_kmers_in_query * prob_unmatched_kmers_not_in_target

    return false_negative_prob 


def compute_bf_false_pos(num_hashes, num_inserted_elements, bf_size):
    '''Computes false positive rate for querying a single element/minimizer'''
    false_pos = (1 - math.exp(
        -1*num_hashes*num_inserted_elements/bf_size
    ))**num_hashes

    return false_pos


def error_lower_bound(error_rate, num_query_minimizers, conf_interval):

    interval_prob = (1.0 - conf_interval)/2

    containment_thresh = error_to_containment(error_rate, kmer_length=16)
    min_hits = max(math.ceil(num_query_minimizers*containment_thresh), 1)

    while min_hits <= num_query_minimizers:
        prob_hit = binom.sf(min_hits-1, num_query_minimizers, 
                            containment_thresh)

        if prob_hit < interval_prob:
            min_hits -= 1   # Last guess was right
            break

        min_hits += 1

    adjusted_containment = min_hits/num_query_minimizers

    lower_bound_error = containment_to_error(adjusted_containment, kmer_length=16)

    return lower_bound_error


def compute_false_hit(false_pos, num_minimizers, error_rate, kmer_length):
    '''Computes the probability of a false bucket hit for a set of minimizers'''
    false_hit_prob = (false_pos)**num_minimizers
    if (error_rate != 0):
        containment = error_to_containment(error_rate, kmer_length)
        num_matches = math.ceil(num_minimizers*containment)
        false_hit_prob = binom.sf(num_matches, num_minimizers, false_pos)
    return false_hit_prob


def compute_sensitivity(true_error_rate, num_query_minimizers, kmer_length):
    true_jaccard_containment = error_to_containment(true_error_rate, 
                                                    kmer_length)

    expected_num_matched = math.ceil(error_to_containment(
        true_error_rate + 0.02, kmer_length
    ) * num_query_minimizers)

    sensitivity = binom.sf(expected_num_matched, num_query_minimizers, 
                           true_jaccard_containment)

    #sensitivity = sensitivity**384

    return sensitivity


def compute_specificity(false_positive, false_negative):
    TN = 1 - false_positive
    return TN / (TN + false_positive)


def compute_bf_size(num_inserted_elements, false_prob_rate):
    '''Computes the array size of the bloom filter given false positive rate
    and number of elements inserted'''
    num_bits = -num_inserted_elements/math.log((1 - false_prob_rate))
    return num_bits


def compute_bf_size_seq(seq_length, window_size, false_prob_rate):
    '''Computes Bloom filter array size from sequence length'''
    num_inserted = compute_num_minimizers(seq_length, window_size)
    bf_size = compute_bf_size(num_inserted, false_prob_rate)

    return bf_size


def compute_bigsi_metrics(bigsi_parameters, query_size, query_sub_rate):
    bf_stats = compute_bigsi_stats(bigsi_parameters)
    bf_false_pos = bf_stats['false positive rate']
    target_size = bigsi_parameters['bin_seq_len']
    window_size = bigsi_parameters['window_size']
    num_query_minimizers = compute_num_minimizers(query_size, window_size)
    #query_sub_rate_lower = error_lower_bound(query_sub_rate, 
    #                                         num_query_minimizers, 
    #                                         conf_interval=0.9999)
    #query_sub_rate = query_sub_rate + (query_sub_rate - query_sub_rate_lower)
    false_pos_prob = compute_false_hit(bf_false_pos, num_query_minimizers, 
                                       query_sub_rate+0.02, 
                                       kmer_length=16)
    false_positive_rate = (1 - (1 - false_pos_prob)**384)

    # false_negative_rate = compute_winnow_false_neg(query_size, target_size, 
    #                                                window_size, query_sub_rate)

    sensitivity = compute_sensitivity(query_sub_rate,  num_query_minimizers, 
                                      kmer_length=16,)
    specificity = 1 - false_positive_rate

    return sensitivity, specificity


def compute_bigsi_stats(parameters):

    kmer_length = parameters['kmer_len']
    max_seq_length = parameters['bin_seq_len']
    window_size = parameters['window_size']
    num_inserted = compute_num_minimizers(max_seq_length, window_size)
    num_minimizers_query = compute_num_minimizers(parameters['min_query_len'],
                                                  window_size)
    error_rate = parameters['error_rate']
    # Take the 90% upper CI of error rate to account for variance in estimate
    # error_rate_lower = error_lower_bound(error_rate, num_minimizers_query, 
    #                                      conf_interval=0.9995)
    # error_rate = error_rate + (error_rate - error_rate_lower)
    # print(error_rate)

    false_hit_thresh = parameters['false_hit_thresh']
    num_cols = parameters['num_cols']

    for bf_size_bits in range(1, 10**8, 10**3):
        bf_size_mb = bits_to_mb(bf_size_bits)
        # num_hashes = compute_num_hash_funcs(bf_size_bits, num_inserted)
        num_hashes = 1
        false_pos_rate = compute_bf_false_pos(num_hashes, num_inserted,
                                              bf_size_bits)
        #j_null = compute_j_null(kmer_length=16, target_length=num_inserted)
        false_hit = compute_false_hit(false_pos_rate,
                                      num_minimizers_query,
                                      error_rate, kmer_length)
        false_hit_total = false_hit*num_cols
        bf_stats = {
            'bf size': bf_size_bits,
            'bigsi size mb': bf_size_mb*num_cols,
            'optimal number hashes': num_hashes,
            'false positive rate': false_pos_rate,
            'false hit rate': false_hit,
            'false hit total': false_hit_total,
        }
        if false_hit_total <= false_hit_thresh:
            return bf_stats


def compute_j_null(kmer_length, target_length):
    return 1 - ((1 - 4**(-1*kmer_length))**target_length)


def compute_prob_false_hit_winnow(j_null, num_query_minimizers, error_rate):

    containment = error_to_containment(error_rate, kmer_length=16)
    num_matches = math.ceil(num_query_minimizers*containment)
    false_hit_prob = binom.sf(num_matches, num_query_minimizers, j_null)

    return false_hit_prob


def compute_minimizer_index_size(seq_length, window_size):
    '''Minimizer index size in mb'''
    minimizer_size = 32 + 16 + 16 + 32  # bits
    num_minimizers = 2*seq_length/window_size
    minimizer_index_size = minimizer_size*num_minimizers/(8*1024*1024)  # mb
    return minimizer_index_size


def compute_containment_diff_bound(delta, window_size, query_size, 
                                   bf_false_pos, true_jaccard_containment):
    scaling_factor = 2/window_size
    num_minimizers_query = scaling_factor*query_size
    exponent_term = (
        (true_jaccard_containment / 
         (true_jaccard_containment + bf_false_pos))**2 
        * delta**2 
        * num_minimizers_query
        * (true_jaccard_containment + bf_false_pos)/3
    )
    return 2*math.exp(-exponent_term)

def main():
    # ar_genes_parameters = {
    #     'bin_seq_len': 7e6,
    #     'window_size': 100,
    #     'min_query_len': 5000,
    #     'error_rate': 0.05,
    #     'false_hit_thresh': 1e-2,
    #     'num_cols': 1000,
    #     'kmer_len': 16,
    # }

    # gene_fusion_parameters = {
    #     'bin_seq_len': 15e6,
    #     'window_size': 50,
    #     'min_query_len': 500,
    #     'error_rate': 0.15,
    #     'false_hit_thresh': 1e-2,
    #     'num_cols': 16*24,
    #     'kmer_len': 16,
    # }

    # human_viruses_parameters = {
    #     'bin_seq_len': 16e6,
    #     'window_size': 50,
    #     'min_query_len': 500,
    #     'error_rate': 0,
    #     'false_hit_thresh': 1e-2,
    #     'num_cols': 16*24,
    #     'kmer_len': 16,
    # }

    # worm = {
    #     'bin_seq_len': 7e6,
    #     'window_size': 100,
    #     'min_query_len': 5000,
    #     'error_rate': 0.05,
    #     'false_hit_thresh': 1e-2,
    #     'num_cols': 16*50,
    #     'kmer_len': 16,
    # }


    # hg38_chr1 = {
    #     'bin_seq_len': 16e6,
    #     'window_size': 100,
    #     'min_query_len': 5000,
    #     'error_rate': 0.15,
    #     'false_hit_thresh': 1e-2,
    #     'num_cols': 16*1,
    #     'kmer_len': 16,
    # }


    # bacterial_ref_sizes = [5682322, 3862530, 4411532, 2821361, 5594605, 4215606, 
    #                    1042519, 2944528, 4042929, 4951383, 4828820, 1641481, 
    #                    4641652, 6264404, 2032807]


    hg38 = {
        'bin_seq_len': 1e6,
        'window_size': 100,
        'min_query_len': 5000,
        'error_rate': 0.05,
        'false_hit_thresh': 1e-2,
        'num_cols': 3102,
        'kmer_len': 16,
    }

    print(compute_minimizer_index_size(3e9, 2000))
    #bigsi_stats = compute_bigsi_stats(hg38)
    #print(hg38, '\n', bigsi_stats)

    # print(compute_containment_diff_bound(0.18, 50, 5000, 0.1749, 0.45))

    # print('hg38')
    # false_negative_prob = compute_winnow_false_neg(
    #     query_size=5000, 
    #     target_size=hg38['bin_seq_len'], 
    #     window_size=hg38['window_size'],
    #     error_rate=0.10
    # )
    # print(false_negative_prob)
    # seq_lengths = [1000, 2000, 3000, 4000, 5000, 10000, 20000, 40000, 80000, 
    #                160000, 200000, 250000, 300000]
    # 
    # query_size_sensitivities = []
    # query_size_specificities = []
    # for query_size in seq_lengths:
    #     sensitivity, specificity = compute_bigsi_metrics(hg38, query_size, 
    #                                     query_sub_rate=0.0)
    #     query_size_sensitivities.append(sensitivity)
    #     query_size_specificities.append(specificity)

    # error_rates = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1]
    # error_rate_sensitivities = []
    # error_rate_specificities = []
    # for query_sub_rate in error_rates:
    #     sensitivity, specificity = compute_bigsi_metrics(hg38, 10000, 
    #                                                      query_sub_rate)
    #     error_rate_sensitivities.append(sensitivity)
    #     error_rate_specificities.append(specificity)

    # print('error_rate_theory_sensitivities =', error_rate_sensitivities)
    # print('error_rate_theory_specificities =', error_rate_specificities)
    # print('query_size_theory_sensitivities =', query_size_sensitivities)
    # print('query_size_theory_specificities =', query_size_specificities)
    # j_null = compute_j_null(kmer_length=16, target_length=1e7)
    # num_query_minimizers = compute_num_minimizers(seq_length=1000, 
    #                                               window_size=100)
    # false_hit_winnow = compute_prob_false_hit_winnow(j_null, 
    #                                                  num_query_minimizers, 
    #                                                  error_rate=0.05)

    # print(j_null, false_hit_winnow)
 
main()
