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


def get_read_by_name(bam: str, read_name: str, loc: dict):
    '''Gets a read filtered by its name and location'''

    with pysam.AlignmentFile(bam, "rb") as samfile:
        for read in samfile.fetch(loc['ref'], loc['start'], loc['end']):
            if read.query_name == read_name:
                return read


def get_sequence(identifier, start, end):
    '''Retrieves sequence from NCBI'''
    Entrez.email = 'shihabdider@berkeley.edu'
    fasta_handle = Entrez.efetch(db='nucleotide', id=identifier, 
                                 rettype='fasta', retmode='text', 
                                 seq_start=start, seq_stop=end)

    fasta_record = SeqIO.read(fasta_handle, "fasta")
    fasta_handle.close()

    return str(fasta_record.seq)


def get_gene_sequence(gene_id):
    '''Retrieves sequence of gene from NCBI'''
    Entrez.email = 'shihabdider@berkeley.edu'
    # get the ref length
    gene_handle = Entrez.efetch(db='gene', id=gene_id,
                           rettype='gene_table', retmode='text')

    #gene_record = SeqIO.read(gene_handle, "text")

    for line in gene_handle:
        print(line)
    gene_handle.close()


def get_random_sequence(identifier, query_len):
    '''Retrieves a random sequence record of specified length from NCBI'''

    Entrez.email = 'shihabdider@berkeley.edu'
    # get the ref length
    gb_handle = Entrez.efetch(db='nucleotide', id=identifier,
                              rettype='gb', retmode='text')
    gb_record = SeqIO.read(gb_handle, "genbank")
    gb_handle.close()
    ref_len = len(gb_record.seq)

    # get random sequence
    query_start = random.randint(0, ref_len - query_len)
    query_end = query_start + query_len - 1
    fasta_handle = Entrez.efetch(db='nucleotide', id=identifier, 
                                 rettype='fasta', retmode='text', 
                                 seq_start=query_start, seq_stop=query_end)

    fasta_record = SeqIO.read(fasta_handle, "fasta")
    fasta_handle.close()

    return fasta_record


def load_query_file(query_path):
    query_record = SeqIO.read(query_path, "fasta")
    return str(query_record.seq)


def run_bigsi_query(query_seq, config):
    '''Runs the bigsi query for a specified bigsi matrix'''

    bigsi = config['bigsi']
    bigsi_path = config['bigsi_path']
    bigsi_config = config['bigsi_config_path']

    query_bigsi_cmd = (
        r"node {0}"
        " -s {1} -b {2} -c {3}").format(
            bigsi, query_seq, bigsi_path, bigsi_config)
    with subprocess.Popen(query_bigsi_cmd,
                          stdout=subprocess.PIPE, shell=True) as proc:
        output = proc.stdout.read().decode('utf-8')
        mappings = output.split('\n')
        return mappings


def get_species_seqs(seq_ids, seq_length, num_queries):
    '''
        Gets random subsequence records from NCBI genome

        Args:
            seq_ids - array of NCBI formatted sequence ids
            seq_length - length of the sequence to be retrieved
            num_queries - how many subsequences to retrieve for each sequence
            id

        Returns:
            records - array of SeqIO records
    '''
    records = []
    for _ in range(num_queries):
        for seq_id in seq_ids:
            random_query_record = get_random_sequence(
                seq_id,
                seq_length
            )
            records.append(random_query_record)

    return records


def run_species_benchmark(benchmark_params):
    '''
    Input:
        benchmark_params = {
            'species_name',
            'seq_ids',
            'query_len',
            'num_queries',
            'bigsi_path',
            'bigsi_config_path',
        }
    Output: runs query_bigsi on benchmark data
    '''

    records = get_species_seqs(
        benchmark_params['seq_ids'],
        benchmark_params['query_len'],
        benchmark_params['num_queries']
    )

    mappings = []
    for record in records:
        query_output = run_bigsi_query(str(record.seq))
        seq_ref = str(record.id).split(':')[0]
        query_start, query_end = str(record.id).split(':')[1].split('-')
        mapping = '\t'.join(
            [
                benchmark_params['species_name'],
                seq_ref,
                query_start,
                query_end,
                str(int(query_end) - int(query_start)),
                query_output.replace('\n', ',')
            ]
        )
        print(mapping)
        mappings.append(mapping)

    return mappings, records


def records_to_fasta(records, output):
    with open(output, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        print('{} reads saved to {}'.format(len(records), output))


def run_mashmap(query, config, identity, output='mashmap.out'):
    '''Runs mashmap on a set of query seqs vs. ref and outputs to file'''

    mashmap_cmd = (
        r"{0}"
        " -q {1} -r {2} -o {3}"
        " -s {4} --pi {5}"
    ).format(config['mashmap'], query, config['ref'], output, 
             config['seq_length'], identity)

    p = subprocess.Popen(mashmap_cmd, shell=True)
    p.communicate()


def get_random_bigsi_bin():
    '''Uses the bigsi bin mapping to retrieve the bounds of a random bin'''

    bigsi_bin_mapping = "../bigsis/hg38_whole_genome_005_bucket_map.json"

    with open(bigsi_bin_mapping, "r") as read_file:
        bin_mapping = json.load(read_file)
        random_bin = random.choice(list(bin_mapping.keys()))
        return random_bin, bin_mapping[random_bin]


def get_bigsi_bin(bin_num: int) -> dict:
    '''
        Uses the bigsi bin mapping file to retrieve the bounds of a given bin

        args:
            bin_num - the bin to retrieve (ranges from 0 to 383)
    '''
    bigsi_bin_mapping = "../bigsis/hg38_whole_genome_005_bucket_map.json"

    with open(bigsi_bin_mapping, "r") as read_file:
        bin_mapping = json.load(read_file)
        return bin_mapping[str(bin_num)]

def run_pacbio_benchmark() -> list:
    '''Runs bigsi queries on Ultralong Nanopore long reads aligned to GRCh38
    reference.'''

    pacbio_longreads = (
        "https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/PacBio_SequelII_CCS_11kb/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam"
    )

    #random_bin = get_bigsi_bin(6)
    reads = []
    reads_per_bin = []
    for i in range(383):
        bigsi_bin = get_bigsi_bin(i)
        gap_width = 10000
        rand_start = random.randint(0, bigsi_bin['bucketEnd'] - gap_width)
        rand_end = rand_start + gap_width - 1

        acn_convert_df = pd.read_csv('./hg38_acn_conversion.txt', sep='\t', 
                                    header=0)
        ref_name = acn_convert_df[
            acn_convert_df['RefSeq-Accn'].str.contains(bigsi_bin['refName'])
        ]['# Sequence-Name'].values[0]

        rand_loc = {
            'ref': ref_name,
            'start': rand_start,
            'end': rand_end,
        }

        aligned_reads = get_aligned_reads(pacbio_longreads, rand_loc)
        reads_per_bin.append(len(aligned_reads))
        reads += aligned_reads

    #print('total num reads: ', len(reads))
    #print('reads per bin: ', reads_per_bin)

    mappings = []
    if len(reads) > 0:
        for read in reads:
            query_output = run_bigsi_query(read.query_alignment_sequence)
            mapping = '\t'.join(
                [read.query_name,
                 read.reference_name, 
                 str(read.reference_start), 
                 str(read.reference_end),
                 str(read.query_alignment_length), 
                 query_output.replace('\n', ',')])
            print(mapping)
            mappings.append(mapping)

        return mappings
    else:
        print('No reads in region', ref_name, rand_start, rand_end)

def run_nanopore_benchmark() -> list:
    '''Runs bigsi queries on Ultralong Nanopore long reads aligned to GRCh38
    reference.'''

    nanopore_longreads = (
        "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"
    )

    #random_bin = get_bigsi_bin(6)
    reads = []
    for i in range(383):
        bigsi_bin = get_bigsi_bin(i)
        gap_width = 10000
        rand_start = random.randint(0, bigsi_bin['bucketEnd'] - gap_width)
        rand_end = rand_start + gap_width - 1

        acn_convert_df = pd.read_csv('./hg38_acn_conversion.txt', sep='\t', 
                                    header=0)
        ref_name = acn_convert_df[
            acn_convert_df['RefSeq-Accn'].str.contains(
                bigsi_bin['refName'])]['UCSC-style-name'].values[0]

        rand_loc = {
            'ref': ref_name,
            'start': rand_start,
            'end': rand_end,
        }

        reads += get_aligned_reads(nanopore_longreads, rand_loc)

    #print('total num reads: ', len(reads))

    mappings = []
    if len(reads) > 0:
        for read in reads:
            query_output = run_bigsi_query(read.query_alignment_sequence)
            mapping = '\t'.join(
                [read.query_name,
                 read.reference_name, 
                 str(read.reference_start), 
                 str(read.reference_end),
                 str(read.query_alignment_length), 
                 query_output.replace('\n', ',')])
            print(mapping)
            mappings.append(mapping)

        return mappings
    else:
        print('No reads in region', ref_name, rand_start, rand_end)


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


def get_multibin_reads(bamfile, mappings):
    '''
        Retrieves the reads which were found in multiple bins after doing the 
        BIGSI query
    '''

    multibin_read_data = []
    for mapping in mappings:
        mapping_list = mapping.split('\t')
        read_info = mapping_list[0:4]
        bigsi_output = mapping_list[-1]
        if bigsi_output:
            bigsi_mappings = bigsi_output.split(',')[0:-1]
            if len(bigsi_mappings) > 1:
                multibin_read_data.append(read_info)

    multibin_reads = []
    for read_info in multibin_read_data:
        read_name = read_info[0]
        read_loc = {
            'ref': read_info[1],
            'start': int(read_info[2]),
            'end': int(read_info[3]),
        }
        read = get_read_by_name(bamfile, read_name, read_loc)
        multibin_reads.append(read)

    return multibin_reads

def reads_to_fasta(reads, output):
    '''Saves reads to FASTA file'''

    records = []
    for read in reads:
        record = SeqRecord(
                    Seq(read.query_alignment_sequence),
                    id=read.query_name,
                    name='{0} {1} {2}'.format(
                        read.reference_name,
                        read.reference_start,
                        read.reference_end),
                )
        records.append(record)

    with open(output, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        print('{} reads saved to {}'.format(len(records), output))


def mappings_to_dict(map_file):
    '''Loads Mashmap output into a dictionary'''
    with open(map_file, 'r') as handle:
        mapping_dict = {}
        for line in handle:
            hit = line.split(' ')
            key = hit[0]
            value = ' '.join([hit[5], hit[7], hit[8]])
            if key not in mapping_dict:
                mapping_dict[key] = [value]
            else:
                mapping_dict[key].append(value)

    return mapping_dict


def simulation_benchmark():
    hg38 = '../../../seqs/human/hg38/ncbi-genomes-2021-11-16/hg38.fna'

    hg38_ids = [
        'NC_000001.11',
        'NC_000002.12',
        'NC_000003.12',
        'NC_000004.12',
        'NC_000005.10',
        'NC_000006.12',
        'NC_000007.14',
        'NC_000008.11',
        'NC_000009.12',
        'NC_000010.11',
        'NC_000011.10',
        'NC_000012.12',
        'NC_000013.11',
        'NC_000014.9',
        'NC_000015.10',
        'NC_000016.10',
        'NC_000017.11',
        'NC_000018.10',
        'NC_000019.10',
        'NC_000020.11',
        'NC_000021.9',
        'NC_000022.11',
        'NC_000023.11',
        'NC_000024.10',
    ]

    hg38_benchmark = {
        'species_name': 'homo sapiens',
        'seq_ids': hg38_ids,
        'query_len': 10000,
        'num_queries': 1,
    }

    hg38_mappings, hg38_records = run_species_benchmark(
        hg38_benchmark
    )

    if hg38_mappings and hg38_records:
        records_to_fasta(hg38_records, 'hg38_random_seqs.fasta')
        run_mashmap('hg38_random_seqs.fasta', hg38, 'human_human.out', 
                    error_rate=100)
        hg38_mapping_dict = mappings_to_dict('human_human.out')
        hg38_sensitivity = compute_sensitivity_species(hg38_mappings,
                                                       hg38_mapping_dict)
        print('human sensitivity: ', hg38_sensitivity)


def gorilla_benchmark():
    hg38 = '../../../seqs/human/hg38/ncbi-genomes-2021-11-16/hg38.fna'

    gorilla_ids = ['NC_044{0}.1'.format(num) for num in range(602, 626)]
    gorilla_benchmark = {
        'species_name': 'gorilla gorilla',
        'seq_ids': gorilla_ids,
        'query_len': 10000,
        'num_queries': 10,
    }

    gorilla_mappings, gorilla_records = run_species_benchmark(
        gorilla_benchmark
    )
    if gorilla_mappings and gorilla_records:
        records_to_fasta(gorilla_records, 'gorilla_random_seqs.fasta')
        run_mashmap('gorilla_random_seqs.fasta', hg38, 'gorilla_human.out')
        gorilla_mapping_dict = mappings_to_dict('gorilla_human.out')
        gorilla_sensitivity = compute_sensitivity_species(gorilla_mappings,
                                                          gorilla_mapping_dict)
        print('gorilla sensitivity: ', gorilla_sensitivity)


def chimp_benchmark():
    hg38 = '../../../seqs/human/hg38/ncbi-genomes-2021-11-16/hg38.fna'

    chimp_ids = ['NC_036{0}.1'.format(num) for num in range(879, 903)]

    chimp_benchmark = {
        'species_name': 'pan trog',
        'seq_ids': chimp_ids,
        'query_len': 10000,
        'num_queries': 10,
    }

    chimp_mappings, chimp_records = run_species_benchmark(chimp_benchmark)
    if chimp_mappings and chimp_records:
        records_to_fasta(chimp_records, 'chimp_random_seqs.fasta')
        run_mashmap('chimp_random_seqs.fasta', hg38, 'chimp_human.out')
        chimp_mapping_dict = mappings_to_dict('chimp_human.out')
        chimp_sensitivity = compute_sensitivity_species(chimp_mappings,
                                                        chimp_mapping_dict)
        print('chimp sensitivity: ', chimp_sensitivity)


def nanopore_benchmark():
    nanopore_mappings = run_nanopore_benchmark()
    if nanopore_mappings:
        nanopore_sensitivity = compute_sensitivity(nanopore_mappings,
                                                   'UCSC-style-name')
        print('nanopore sensitivity: ', nanopore_sensitivity)
        #nanopore_multibin_reads = get_multibin_reads(nanopore_longreads,
        #                                             nanopore_mappings)
        #reads_to_fasta(nanopore_multibin_reads,
        #               'multibin_reads_nanopore.fasta')


def pacbio_benchmark():
    pacbio_mappings = run_pacbio_benchmark()
    if pacbio_mappings:
        pacbio_sensitivity = compute_sensitivity(pacbio_mappings,
                                                 '# Sequence-Name')
        print('pacbio sensitivity: ', pacbio_sensitivity)
        #pacbio_multibin_reads = get_multibin_reads(pacbio_longreads,
        #                                           pacbio_mappings)
        #reads_to_fasta(pacbio_multibin_reads,
        #               'multibin_reads_pacbio.fasta')


def write_to_json(data, output):
    with open(output, 'w') as handle:
        json_data = json.dumps(data)
        handle.write(json_data)

def main():
    parser = argparse.ArgumentParser(
        description="Run Flashmap on benchmark sequences.")
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
        "-i", "--identity", type=int, 
        help="Identity cutoff for mashmap", 
        required=True
    )
    parser.add_argument(
        "-o", "--output", type=str, 
        help="Base file path for outputting benchmark files", 
        required=True
    )
    args = parser.parse_args()

    config = ''
    with open(args.config, 'r') as handle:
        config = json.load(handle)

    query_records = [record for record in SeqIO.parse(args.query, 'fasta')]

    bigsi_results = {}
    for record in query_records:
        query_seq = str(record.seq)
        bigsi_output = run_bigsi_query(query_seq, config)
        bigsi_results[record.id] = bigsi_output

    bigsi_results_path = args.output + '.bigsi.json'
    write_to_json(bigsi_results, bigsi_results_path)

    mashmap_results_path = args.output + '.mashmap.out'
    run_mashmap(args.query, config, args.identity, output=mashmap_results_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
