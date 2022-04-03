from Bio import Entrez
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import random
import subprocess
import pysam
import json
import pandas as pd




def get_aligned_reads(
    bam: str, loc: dict, identity_threshold: int = 0.95
) -> list:
    '''
        Fetches read sequences from a bam alignment (local or remote) that
        align to reference at loc with identity higher than threshold.

        args:
            bam - path to local or remote (e.g ftp) bam/sam/cram file
            ref
            loc - location in reference within which read aligns, consists of 
            ref_name, start and end
            identity_threshold - BLAST-like identity score between read and 
            reference

        returns:
            reads - array of read objects

    '''

    reads = []
    with pysam.AlignmentFile(bam, "rb") as samfile:
        for read in samfile.fetch(loc['ref'], loc['start'], loc['end']):
            is_query_right_size = (read.query_alignment_length > 5000 and
                                   read.query_alignment_length < 300000)
            is_mapped = not read.is_unmapped
            is_good_quality = (read.mapping_quality >= 20 and 
                               read.mapping_quality != 255)
            num_matches = len(read.get_aligned_pairs(matches_only=True))
            total_seq = len(read.get_aligned_pairs())
            is_identity = (num_matches/total_seq) > identity_threshold
            if (read.query_alignment_sequence 
                    and is_query_right_size and is_mapped and is_good_quality 
                    and is_identity):
                reads.append(read)

    return reads


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
    '''Retrieves a random sequence of specified length from NCBI'''

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

    return str(fasta_record.seq)


def load_query_file(query_path):
    query_record = SeqIO.read(query_path, "fasta")
    return str(query_record.seq)


def run_bigsi_query(query_seq):
    '''Runs the bigsi query for a specified bigsi matrix'''

    bigsi_path = '../bigsis/hg38_whole_genome_005.bin'
    bigsi_config_path = '../bigsi.random.query.config.json'

    query_bigsi_cmd = (
        r"node ../bin/query_bigsi.js"
        " -s {0} -b {1} -c {2}").format(
            query_seq, bigsi_path, bigsi_config_path)
    with subprocess.Popen(query_bigsi_cmd,
                          stdout=subprocess.PIPE, shell=True) as proc:
        output = proc.stdout.read().decode('utf-8')
        return output


def run_species_benchmark(benchmark_params):
    '''
    Input:
        benchmark_params = {
            'species_name',
            'record_id',
            'query_len',
            'num_queries',
            'bigsi_path',
            'bigsi_config_path',
        }
    Output: runs query_bigsi on benchmark data
    '''

    for i in range(benchmark_params['num_queries']):
        random_query_seq = get_random_sequence(benchmark_params['record_id'],
                                               benchmark_params['query_len'])
        run_bigsi_query(random_query_seq, 
                        benchmark_params['bigsi_path'], 
                        benchmark_params['bigsi_config_path'])


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

def run_nanopore_benchmark() -> list:
    '''Runs bigsi queries on Ultralong Nanopore long reads aligned to GRCh38
    reference.'''

    nanopore_longreads = (
        #"https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/PacBio_SequelII_CCS_11kb/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam"
        "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"
    )

    #random_bin = get_bigsi_bin(6)
    reads = []
    for i in range(16):
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

    print('total num reads: ', len(reads))

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


def compute_sensitivity(mappings):
    '''Given a list of mappings computes the sensitivity of bigsi search'''
    true_positives = 0
    total = 0
    acn_convert_df = pd.read_csv('./hg38_acn_conversion.txt', sep='\t', 
                                 header=0)
    num_buckets_hit = []
    for mapping in mappings:
        is_true_positive = False
        mapping_list = mapping.split('\t')
        read_ref = mapping_list[0]
        read_start = int(mapping_list[1])
        read_end = int(mapping_list[2])
        bigsi_output = mapping_list[-1]
        if bigsi_output:
            bigsi_mappings = bigsi_output.split(',')[0:-1]
            num_buckets_hit.append(len(bigsi_mappings))
            for bigsi_map in bigsi_mappings:
                bigsi_map_list = bigsi_map.split(' ')
                mapped_ref_name = acn_convert_df[
                    acn_convert_df['RefSeq-Accn'].str.contains(
                        bigsi_map_list[0])]['UCSC-style-name'].values[0]
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
    print(num_buckets_hit)
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


def main():
    #random.seed(1)

    chimp_benchmark = {
        'species_name': 'pan trog',
        'record_id': 'NC_036884.1',
        'query_len': 300000,
        'num_queries': 20,
    }


    #get_gene_sequence('454734')
    #run_species_benchmark(chimp_benchmark)
    #random_sequence = get_random_sequence('NC_036883.1', 300000)
    #chimp_gene = get_sequence('NC_036896.1', 10480693, 10531554)
    #dog_gene = get_sequence('NC_051813.1', 19428782, 19464638)
    #gene = get_sequence('NC_000086.8', 162922338, 162971414)
    #run_bigsi_query(gene, bigsi_path, bigsi_config_path)
    nanopore_longreads = (
        #"https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/NA12878/PacBio_SequelII_CCS_11kb/HG001.SequelII.pbmm2.hs37d5.whatshap.haplotag.RTG.trio.bam"
        "ftp://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/"
        "NA12878/Ultralong_OxfordNanopore/NA12878-minion-ul_GRCh38.bam"
    )
    mappings = run_nanopore_benchmark()
    if mappings:
        multibin_reads = get_multibin_reads(nanopore_longreads, mappings)
        reads_to_fasta(multibin_reads, 'multibin_reads.fasta')
    #print(get_random_bigsi_bin('../bigsis/hg38_whole_genome_005_bucket_map.json'))
    #acn_convert_df = pd.read_csv('./hg38_acn_conversion.txt', sep='\t', 
    #                             header=0)
    #contain_values = acn_convert_df[acn_convert_df['UCSC-style-name'].str.contains('chr1')]
    #print(contain_values)

main()
