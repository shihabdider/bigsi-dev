'''
Complexity calculations for BIGSI and Mashmap

Questions:
    - What are all the parameters involved in constructing the sketch bigsi?
        - Reference size
        - Containment threshold
        - Sensitivity threshold
        - Bloom filter length
        - Number of elements inserted into Bloom filter
        - Number of bins
        - Length of sequence for each bin
        - Window size for sketching
        - Min fragment size for fragment search
    - Which of those parameters are exposed?
        - Reference size
        - Containment threshold
        - Sensitivity threshold
        - Window size
        - Min fragment size
        - Number of bins
    - What is the relationship between the internal parameters and those
      exposed to the user?
    - How do you compute the sensitivity of a query?
    - How do you interpret a query? In terms of the size of the sequences?
    - What are some possible irl queries?
        - Anti-biotic resistance genes
        - Viruses in human genomes
        - Gene fusions
        - Segmental duplications
        - Cancer copy number variants
        - HIV integration sites
    - What are the sizes of possible irl queries?
        - For antibiotic resistance genes:
            - Source: CARD
            - Range: 162-4359bp
            - Average (std): 957 (320)
            - Number of genomes: 586
            - Number of representative bacterial sequences in RefSeq: 14,906
            - Bigsi sizes (linear in both resolution and number of genomes; the 
            following assumes fixed 30mb bin sizes and exact matches):
                - 100bp resolution, 1000 genomes: 101Mb
                - 500bp resolution, 1000 genomes: 20Mb
                - 100bp resolution, 15k genomes: 1.5Gb
                - 500bp resolution, 15k genomes: 300Mb
        - For viruses:
            - Source: NCBI
            - Range: 162-4359bp
            - Average (std): 957 (320)
            - Number of genomes: 586
            - Number of representative bacterial sequences in RefSeq: 14,906
            - Bigsi sizes (linear in both resolution and number of genomes; the 
            following assumes fixed 30mb bin sizes and exact matches):
                - 100bp resolution, 1000 genomes: 101Mb
                - 500bp resolution, 1000 genomes: 20Mb
                - 100bp resolution, 15k genomes: 1.5Gb
                - 500bp resolution, 15k genomes: 300Mb
        - For gene fusions:
            - Source: FusionGDB
            - Range: 162-4359bp
            - Average (std): 957 (320)
            - Number of genomes: 586
            - Number of representative bacterial sequences in RefSeq: 14,906
            - Bigsi sizes (linear in both resolution and number of genomes; the 
            following assumes fixed 30mb bin sizes and exact matches):
                - 100bp resolution, 1000 genomes: 101Mb
                - 500bp resolution, 1000 genomes: 20Mb
                - 100bp resolution, 15k genomes: 1.5Gb
                - 500bp resolution, 15k genomes: 300Mb
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


def compute_false_hit(false_pos, num_minimizers, error_rate, kmer_length):
    '''Computes the probability of a false bucket hit for a set of minimizers'''
    false_hit_prob = (false_pos)**num_minimizers
    if (error_rate != 0):
        containment = error_to_containment(error_rate, kmer_length)
        num_matches = math.ceil(num_minimizers*containment)
        false_hit_prob = binom.sf(num_matches, num_minimizers, false_pos)
    return false_hit_prob


def compute_sensitivity(false_positive, false_negative):
    TP = 1 - false_positive
    return TP / (TP + false_negative)


def compute_specificity(false_positive, false_negative):
    TN = 1 - false_negative
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
    print(bf_stats['bf size'])
    false_pos_rate = bf_stats['false positive rate']
    num_query_minimizers = compute_num_minimizers(query_size, window_size=100)
    false_pos_prob = compute_false_hit(false_pos_rate, num_query_minimizers, 
                                       query_sub_rate, kmer_length=16)
    false_pos_prob = (1 - (1 - false_pos_prob)**384)

    target_size = bigsi_parameters['bin_seq_len']
    window_size = bigsi_parameters['window_size']
    false_negative_prob = compute_winnow_false_neg(query_size, target_size, 
                                                   window_size, query_sub_rate)

    sensitivity = compute_sensitivity(false_pos_prob, false_negative_prob)
    specificity = compute_specificity(false_pos_prob, false_negative_prob)

    return sensitivity, specificity


def compute_bigsi_stats(parameters):

    kmer_length = parameters['kmer_len']
    max_seq_length = parameters['bin_seq_len']
    window_size = parameters['window_size']
    num_inserted = compute_num_minimizers(max_seq_length, window_size)
    num_minimizers_query = compute_num_minimizers(parameters['min_query_len'],
                                                  window_size)
    error_rate = parameters['error_rate']
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



ar_genes_parameters = {
    'bin_seq_len': 7e6,
    'window_size': 100,
    'min_query_len': 5000,
    'error_rate': 0.05,
    'false_hit_thresh': 1e-2,
    'num_cols': 1000,
    'kmer_len': 16,
}

gene_fusion_parameters = {
    'bin_seq_len': 15e6,
    'window_size': 50,
    'min_query_len': 500,
    'error_rate': 0.15,
    'false_hit_thresh': 1e-2,
    'num_cols': 16*24,
    'kmer_len': 16,
}

human_viruses_parameters = {
    'bin_seq_len': 16e6,
    'window_size': 50,
    'min_query_len': 500,
    'error_rate': 0,
    'false_hit_thresh': 1e-2,
    'num_cols': 16*24,
    'kmer_len': 16,
}

worm = {
    'bin_seq_len': 7e6,
    'window_size': 100,
    'min_query_len': 5000,
    'error_rate': 0.05,
    'false_hit_thresh': 1e-2,
    'num_cols': 16*50,
    'kmer_len': 16,
}


hg38_chr1 = {
    'bin_seq_len': 16e6,
    'window_size': 100,
    'min_query_len': 5000,
    'error_rate': 0.15,
    'false_hit_thresh': 1e-2,
    'num_cols': 16*1,
    'kmer_len': 16,
}


bacterial_ref_sizes = [5682322, 3862530, 4411532, 2821361, 5594605, 4215606, 
                       1042519, 2944528, 4042929, 4951383, 4828820, 1641481, 
                       4641652, 6264404, 2032807]


def compute_minimizer_index_size(seq_length, window_size):
    '''Minimizer index size in mb'''
    minimizer_size = 32 + 16 + 16 + 32  # bits
    num_minimizers = 2*seq_length/window_size
    minimizer_index_size = minimizer_size*num_minimizers/(8*1024*1024)  # mb
    return minimizer_index_size



# Minimizers

ALPHABET_SIZE = 4
KMER_SIZE = 15
MIN_MATCH_LENGTH = 5*10**3
REF_LENGTH = 3*10**9


def compute_prob_random_kmer_in_set(kmer_length, kmer_set_size):
    '''Computes the probability of a random k-mer appearing in a k-mer set at
    least once, assuming the k-mers are independent'''
    alphabet_size = 4   # for DNA
    random_kmer_match = 1 - \
        (1 - alphabet_size**(-kmer_length))**(kmer_set_size)

    return random_kmer_match


def compute_null_jaccard():
    '''The null Jaccard score is the expected similiarity between two random
    k-mer sets X and Y given one match. Here it is used as the probability of
    an element randomly selected from the union, also being found in the
    intersection'''
    
    # We assume the seq length >> kmer length and equal size ref and query seqs
    random_kmer_match_ref = compute_prob_random_kmer_in_set(
        minimizer_global_parameters['MIN_MATCH_LENGTH']
    )
    random_kmer_match_query = compute_prob_random_kmer_in_set(
        minimizer_global_parameters['MIN_MATCH_LENGTH']
    )

    prob_random_kmer_in_intersection = random_kmer_match_ref * \
        random_kmer_match_query
    prob_random_kmer_in_union = random_kmer_match_ref + \
        random_kmer_match_query - prob_random_kmer_in_intersection

    j_null = prob_random_kmer_in_intersection/prob_random_kmer_in_union
    return j_null


def mash_dist_to_jaccard(mash_distance):
    jaccard = 1.0 / (2.0 * math.exp(KMER_SIZE*mash_distance) - 1.0)
    return jaccard


def jaccard_to_mash_dist(jaccard_score):
    if jaccard_score == 0:
        return 1.0

    if jaccard_score == 1:
        return 0.0

    mash_distance = (-1.0 / KMER_SIZE) * math.log(
        2.0 * jaccard_score/(1+jaccard_score)
    )

    return mash_distance


def estimate_min_hits(perc_identity, sketch_size):
    mash_distance = 1 - perc_identity/100
    jaccard_score = mash_dist_to_jaccard(mash_distance)
    minimum_shared_minimizers = math.ceil(1.0 * sketch_size * jaccard_score)

    return minimum_shared_minimizers


def mash_dist_lower_bound(mash_distance, sketch_size, conf_interval):
    interval_prob = (1.0 - conf_interval)/2

    jaccard_score = mash_dist_to_jaccard(mash_distance)
    min_hits = max(math.ceil(sketch_size*jaccard_score), 1)

    while min_hits <= sketch_size:
        prob_hit = binom.sf(min_hits-1, sketch_size, jaccard_score)

        if prob_hit < interval_prob:
            min_hits -= 1   # Last guess was right

        min_hits += 1

    jaccard = min_hits/sketch_size

    lower_bound = jaccard_to_mash_dist(jaccard)

    return lower_bound

def estimate_min_hits_relaxed(sketch_size, perc_identity):
    search_upper_bound = estimate_min_hits(sketch_size, perc_identity)

    minimum_hits_relaxed = 0
    for shared_hits in range(search_upper_bound, -1, -1):
        jaccard_score = 1.0*shared_hits/sketch_size
        mash_distance = jaccard_to_mash_dist(jaccard_score)

        distance_lower_bound = mash_distance_lower_bound(mash_distance, 
                                                         sketch_size,
                                                         conf_interval)

        distance_upper_bound = 100*(1 - distance_lower_bound)

        if (distance_upper_bound >= perc_identity):
            minimum_hits_relaxed = shared_hits
        else:
            break

        return minimum_hits_relaxed


def estimate_pval(sketch_size, min_identity):
    jaccard_null = compute_null_jaccard()
    min_hits = estimate_min_hits(min_identity, sketch_size)

    pval = None
    if min_hits == 0:
        cdf_complement = 1.0
    else:
        cdf_complement = binom.sf(min_hits-1, sketch_size, jaccard_null)
        pval = minimizer_global_parameters['REF_LENGTH'] * cdf_complement

    return pval


def compute_optimal_window_size(min_identity, min_pval):
    potential_sketch_sizes = [1, 2, 5] 
    potential_sketch_sizes += list(
        range(10, minimizer_global_parameters['MIN_MATCH_LENGTH'], 10)
    )

    optimal_sketch_size = 0
    for sketch_size in potential_sketch_sizes:
        pval = estimate_pval(sketch_size, min_identity)
        if pval <= min_pval:
            optimal_sketch_size = sketch_size
            break

    window_size = 2.0*minimizer_global_parameters['MIN_MATCH_LENGTH'] \
        / optimal_sketch_size

    return window_size


def compute_expected_sketch_size(window_length, num_kmers):
    sketch_size = 2*num_kmers/window_length
    return sketch_size



def test_estimates_random_kmer():
    kmerSpace = math.pow(4,15)
    lengthQuery = 5000
    pX = compute_prob_random_kmer_in_set(lengthQuery, 15)
    pY = 1. / (1. + kmerSpace / lengthQuery);
    pZ = 1 - math.exp(-lengthQuery/kmerSpace)
    pA = lengthQuery/kmerSpace

    print(pX, pY, pZ, pA)
    print('real - jain approx', pX-pY)
    print('real - exp approx', pX-pZ)
    print('real - binom approx', pX-pA)

def containment_false_pos():
    ref_sketch_size = int(1e4)
    query_sketch_size = int(1e2)
    error_rate = 0.2
    kmer_length = 8
    containment_ratio = math.exp(-kmer_length*error_rate)
    min_num_matches = math.ceil(containment_ratio*query_sketch_size)
    print(min_num_matches)

    r = compute_prob_random_kmer_in_set(kmer_length, ref_sketch_size)
    false_pos_rate = binom.sf(min_num_matches, query_sketch_size, r)
    print(false_pos_rate)


def main():
    hg38 = {
        'bin_seq_len': 16e6,
        'window_size': 100,
        'min_query_len': 5000,
        'error_rate': 0.07,
        'false_hit_thresh': 1e-2,
        'num_cols': 16*24,
        'kmer_len': 16,
    }

    # print('hg38')
    # false_negative_prob = compute_winnow_false_neg(
    #     query_size=5000, 
    #     target_size=hg38['bin_seq_len'], 
    #     window_size=hg38['window_size'],
    #     error_rate=0.10
    # )
    # print(false_negative_prob)
    metrics = compute_bigsi_metrics(hg38, query_size=10000, query_sub_rate=0.03)
    print(metrics)

    # j_null = compute_j_null(kmer_length=16, target_length=1e7)
    # num_query_minimizers = compute_num_minimizers(seq_length=1000, 
    #                                               window_size=100)
    # false_hit_winnow = compute_prob_false_hit_winnow(j_null, 
    #                                                  num_query_minimizers, 
    #                                                  error_rate=0.05)

    # print(j_null, false_hit_winnow)
 
main()
