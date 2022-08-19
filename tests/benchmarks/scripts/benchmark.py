from Bio import SeqIO

import subprocess
import json
import argparse
import logging


def run_bigsi_query(query_seq, config, subrate):
    '''Runs the bigsi query for a specified bigsi matrix'''

    bigsi = config['bigsi']
    bigsi_path = config['bigsi_path']
    bigsi_config = config['bigsi_config_path']

    query_bigsi_cmd = (
        r"node {0}"
        " -s {1} -b {2} -c {3} -e {4}").format(
            bigsi, query_seq, bigsi_path, bigsi_config, subrate)
    with subprocess.Popen(query_bigsi_cmd,
                          stdout=subprocess.PIPE, shell=True) as proc:
        output = proc.stdout.read().decode('utf-8')
        #print(output)
        mappings = output.split('\n')
        return mappings


def run_sketched_mashmap(query, config, identity, seq_length, output='mashmap.out'):
    '''Runs mashmap on a set of query seqs vs. sketched ref and outputs to file'''

    mashmap_cmd = (
        r"{0}"
        " -q {1} --rs {2} -o {3}"
        " -s {4} --pi {5} -m 0"
    ).format(config['mashmap_sketch'], query, config['ref_sketch'], output,
             seq_length, identity)

    p = subprocess.Popen(mashmap_cmd, shell=True)
    p.communicate()

def run_mashmap(query, config, identity, seq_length, output='mashmap.out'):
    '''Runs mashmap on a set of query seqs vs. ref and outputs to file'''

    mashmap_cmd = (
        r"{0}"
        " -q {1} -r {2} -o {3}"
        " -s {4} --pi {5}"
    ).format(config['mashmap'], query, config['ref'], output,
             seq_length, identity)

    p = subprocess.Popen(mashmap_cmd, shell=True)
    p.communicate()


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
    parser.add_argument(
        "-m", "--mashmap", type=int,
        help="Flag (0 or 1) for whether to run mashmap on the benchmark data",
        required=True,
        default=1
    )
    args = parser.parse_args()

    config = {}
    with open(args.config, 'r') as handle:
        config = json.load(handle)

    query_records = [record for record in SeqIO.parse(args.query, 'fasta')]

    bigsi_results = {}
    query_lengths = []
    subrate = int(args.identity)/100
    for record in query_records:
        query_seq = str(record.seq)
        query_lengths.append(len(query_seq))
        bigsi_output = run_bigsi_query(query_seq, config, subrate)
        bigsi_results[record.id] = bigsi_output

    bigsi_results_path = args.output + '.bigsi.json'
    write_to_json(bigsi_results, bigsi_results_path)
    print('Wrote to {}'.format(bigsi_results_path))

    if args.mashmap:
        mashmap_results_path = args.output + '.mashmap.out'
        mashmap_query_length = min(query_lengths) - 1
        mashmap_identity = (1 - subrate) * 100
        run_sketched_mashmap(args.query, config, mashmap_identity, mashmap_query_length, output=mashmap_results_path)


if __name__ == "__main__":
    try:
        main()
    except Exception as e:
        logging.exception(e)
