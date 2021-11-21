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
    - What is the relationship between the internal parameters and the those 
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
        - Anti-biotic resistance genes: {'average': 957.2272230718811,
        'minimum': 162, 'maximum': 4359, 'std': 319.67748306799336}
'''

import math
from scipy.stats import binom

# Sketch BIGSI


def compute_num_hash_funcs(bf_size, num_inserted):
    '''Computes optimal number of hash functions for a bloom filter'''
    return math.ceil((bf_size/num_inserted)*math.log(2))


def compute_bf_false_pos(num_hashes, num_inserted_elements, bf_size):
    '''Computes false positive rate for a single element/minimizer'''
    false_pos = (1 - math.exp(
        -1*num_hashes*num_inserted_elements/bf_size
    ))**num_hashes
    return false_pos


def compute_false_hit(false_pos, num_minimizers, perc_identity):
    '''Computes the probability of a false bucket hit for a query'''
    false_hit_prob = false_pos**num_minimizers
    if (perc_identity != 1.0):
        num_matches = num_minimizers*perc_identity
        false_hit_prob = binom.sf(num_matches, num_minimizers, false_pos)
    return false_hit_prob


def compute_bf_size(false_prob_rate, num_inserted_elements):
    '''Computes the array size of the bloom filter given false positive rate
    and number of elements inserted'''
    num_bits = -num_inserted_elements/math.log((1 - false_prob_rate))
    return num_bits

def bits_to_mb(bits):
    return bits/(8*10**6)

def compute_num_minimizers(seq_length, window_size):
    num_minimizers = seq_length/window_size * 2
    return num_minimizers

def print_bf_stats():
    max_seq_length = 30e6
    window_size = 13
    num_inserted = compute_num_minimizers(max_seq_length, window_size)
    num_minimizers_query = 15
    perc_identity = 1.0
    false_hit_thresh = 1e-2
    num_cols = 100
    for bf_size_bits in range(1, 10**9, 10**3):
        bf_size_mb = bits_to_mb(bf_size_bits)
        num_hashes = compute_num_hash_funcs(bf_size_bits, num_inserted)
        false_pos_rate = compute_bf_false_pos(num_hashes, num_inserted,
                                              bf_size_bits)
        false_hit = compute_false_hit(false_pos_rate,
                                      num_minimizers_query,
                                      perc_identity)
        false_hit_total = false_hit*1
        bf_stats = 'bf size (in bits): {0} \n \
                \t(in Mb) {1} \n \
                optimal number hashes {2} \n \
                false positive rate: {3} \n \
                false hit rate: {4}\n \
                false hit total: {5}'.format(
                bf_size_bits, bf_size_mb*num_cols, num_hashes, 
                false_pos_rate, false_hit, false_hit_total)
        if false_hit_total <= false_hit_thresh:
            print('optimal params found!')
            print(bf_stats)
            break


parameters = {
    'ref_len': 3e8,
    'containment_thresh': 1.0,
    'senstivity_thresh': 0.99,
    'window_size': 0,
    'min_fragment_len': 100,
    'num_bins': 128,
    'kmer_len': 16,
}




def print_bigsi_benchmark(parameters):
    # frag query
    false_pos_rate = 0.4
    num_minimizers = 14
    perc_identity = 1.0
    false_hit_rate = compute_false_hit(false_pos_rate, num_minimizers, perc_identity)
    num_seqs = 24
    num_bigsi_cols = 128*num_seqs
    bigsi_false_hit_rate = 1 - (1 - false_hit_rate)**num_bigsi_cols

    print(false_hit_rate, bigsi_false_hit_rate)


def main():
    print_bf_stats()
    #print_bigsi_benchmark(parameters)
 
main()

# Minimizers

ALPHABET_SIZE = 4
KMER_SIZE = 15
MIN_MATCH_LENGTH = 5*10**3
REF_LENGTH = 3*10**9

def compute_prob_random_kmer_in_set(seq_length):
    '''Computes the probability of a random k-mer appearing in a k-mer set at
    least once, assuming the k-mers are independent (i.e the apperance of one 
    k-mer does not bias the appearence of another)'''
    alphabet_size = 4   # for DNA
    random_kmer_match = 1 - (1 - alphabet_size**(-KMER_SIZE))**(seq_length)
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
