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


def get_read_by_name(bam, read_name, loct):
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



def get_aligned_reads_records(bam, loc, identity_threshold = 0.95):
    '''
        Fetches read sequences from a bam alignment (local or remote) that
        align to reference at loc with identity higher than threshold.

        args:
            bam - path to local or remote (e.g ftp) bam/sam/cram file
            ref
            loc - location in reference within which read aligns, consists of
            ref, start and end
            identity_threshold - BLAST-like identity score between read and
            reference

        returns:
            reads - array of SeqRecords

    '''

    reads_records = []
    with pysam.AlignmentFile(bam, "rb") as bamfile:
        for read in bamfile.fetch(loc['ref'], loc['start'], loc['end']):
            is_right_size_query = (read.query_alignment_length > 5000 and
                                   read.query_alignment_length < 300000)
            is_mapped = not read.is_unmapped
            is_good_quality = (read.mapping_quality >= 20 and
                               read.mapping_quality != 255)
            num_matches = len(read.get_aligned_pairs(matches_only=True))
            total_seq = len(read.get_aligned_pairs())
            is_above_identity_thresh = (num_matches/total_seq) > identity_threshold
            if (read.query_alignment_sequence
                and is_right_size_query
                and is_mapped
                and is_good_quality
                and is_above_identity_thresh):
                read_record = SeqRecord(Seq(read.query_sequence),
                                        read.query_name)
                reads_records.append(read_record)

    return reads_records


def get_example_ref_name(bam):
    '''Retrieves an example ref name from a bam file to determine accession
    style'''

    with pysam.AlignmentFile(bam, "rb") as bamfile:
        # get only the first read
        for read in bamfile.fetch(until_eof=True):
            return read.reference_name


def convert_acn(identifier, target_ref_name_example, acn_file):
    '''Converts identifier to accession style of a given example'''

    # load the table of acns
    acn_convert_df = pd.read_csv(acn_file, sep='\t', header=0)

    found = acn_convert_df.isin([identifier]).any()
    identifier_column = found[found].index.values[0]

    found = acn_convert_df.isin([target_ref_name_example]).any()
    target_ref_name_column = found[found].index.values[0]

    # convert identifier to the target ref name style using the above column
    ref_name = acn_convert_df[
        acn_convert_df[identifier_column].str.contains(identifier)
        ][target_ref_name_column].values[0]

    return ref_name

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


def get_pysam_record(fasta_path, identifier, start, end):
    ref = pysam.FastaFile(fasta_path)
    seq = ref.fetch(identifier, start, end).upper()
    record_id = '{0}:{1}-{2}'.format(identifier, start, end)
    fasta_record = SeqRecord(Seq(seq), record_id)
    return fasta_record


def get_fasta_record(fasta_handle, identifier, start, end):
    for record in SeqIO.parse(fasta_handle, 'fasta'):
        if identifier == record.id:
            sequence = str(record.seq[start:end])
            fasta_record = SeqRecord(Seq(sequence), record.id)
            return fasta_record


def get_ncbi_record(identifier, start, end):
    '''Retrieves seq record from NCBI'''
    Entrez.email = 'shihabdider@berkeley.edu'
    fasta_handle = Entrez.efetch(db='nucleotide', id=identifier, 
                                 rettype='fasta', retmode='text', 
                                 seq_start=start, seq_stop=end)

    fasta_record = SeqIO.read(fasta_handle, "fasta")
    fasta_handle.close()

    return fasta_record


def make_random_fasta_interval(faidx, identifier, interval_length):
    ref_lengths = {}
    with open(faidx, 'r') as handle:
        for line in handle:
            split_line = line.split('\t')
            name = split_line[0]
            length = int(split_line[1])
            ref_lengths[name] = length

    ref_length = ref_lengths[identifier]
    # get the interval
    interval_start = random.randint(0, ref_length - interval_length)
    interval_end = interval_start + interval_length - 1

    return interval_start, interval_end


def make_random_interval(identifier, interval_length):
    '''Generates random start-end for retriving from NCBI'''

    # get the ref length
    Entrez.email = 'shihabdider@berkeley.edu'
    gb_handle = Entrez.efetch(db='nucleotide', id=identifier,
                              rettype='gb', retmode='text')
    gb_record = SeqIO.read(gb_handle, "genbank")
    gb_handle.close()
    ref_length = len(gb_record.seq)

    # get the interval
    interval_start = random.randint(0, ref_length - interval_length)
    interval_end = interval_start + interval_length - 1

    return interval_start, interval_end


def records_to_fasta(records, output):
    with open(output, 'w') as handle:
        SeqIO.write(records, handle, 'fasta')
        print('{} seq records saved to {}'.format(len(records), output))


def main():
    parser = argparse.ArgumentParser(
        description="Create sequence FASTA for running benchmarks"
    )
    parser.add_argument(
        "-a", "--acn", type=str,
        help="path to text file with accession ids",
        required=False
    )
    parser.add_argument(
        "-b", "--bam", type=str, 
        help="path or url to bam alignment file", 
        required=False
    )
    parser.add_argument(
        "-f", "--fasta", type=str,
        help="path to fasta file",
        required=False
    )
    parser.add_argument(
        "-x", "--faidx", type=str, 
        help="path reference faidx (to get ref lengths)", 
        required=True
    )
    parser.add_argument(
        "-i", "--identifiers", type=str, 
        help="path to file containing NCBI identifiers (one per line)", 
        required=True
    )
    parser.add_argument(
        "-l", "--length", type=float,
        help=("length of the sequences"),
        default=10000,
    )
    parser.add_argument(
        "-n", "--num", type=int,
        help="number of sequences to generate per identifier",
        default=1
    )
    parser.add_argument(
        "-o", "--output", type=str,
        help="output path to save FASTA file",
        required=True
    )
    args = parser.parse_args()

    identifiers = []
    with open(args.identifiers, 'r') as handle:
        for line in handle:
            split_line = line.split('\t')
            identifier = split_line[0].rstrip()
            identifiers.append(identifier)

    output_records = []
    if args.fasta:
        print('Retrieving random records from FASTA...')
        with open(args.fasta, 'r') as handle:
            for identifier in identifiers:
                for i in range(args.num):
                    rand_start, rand_end = make_random_fasta_interval(
                        args.faidx, identifier, args.length)
                    random_record = get_pysam_record(args.fasta, identifier,
                                                    rand_start, rand_end)
                    num_n = str(random_record.seq).count('N')
                    is_valid_seq = num_n/len(random_record.seq) < 0.3
                    #print(identifier, i, random_record)
                    if is_valid_seq:
                        output_records.append(random_record)

    elif args.bam and args.acn:
        bam_ref_name_example = get_example_ref_name(args.bam)
        for identifier in identifiers:
            for i in range(args.num):
                rand_start, rand_end = make_random_fasta_interval(
                    args.faidx, identifier, args.length)

                ref_name = convert_acn(identifier, bam_ref_name_example, args.acn)

                rand_loc = {
                    'ref': ref_name,
                    'start': rand_start,
                    'end': rand_end
                }

                output_records += get_aligned_reads_records(args.bam, rand_loc)
    else:
        print('Please provide a FASTA or BAM file')

    records_to_fasta(output_records, args.output)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
