#!/usr/bin/env python3

import argparse
import sys

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def parse_args(args=None):
  Description='Count number of samples for IQ-Tree'

  parser = argparse.ArgumentParser(description=Description)
  parser.add_argument('compiled_fasta_file',
                      help='Complete fasta file that needs to be counted.')
  return parser.parse_args(args)

def count_samples(compiled_fasta_file):
    with open(compiled_fasta_file, "r") as inFasta, open("count.txt", "w") as count_file:
        # Starts at -1 to account for reference fasta 
        count = -1
        for each_record in SeqIO.parse(inFasta, "fasta"):
            count += 1
        count_file.write(str(count))
    return count_file

def main(args=None):
    args = parse_args(args)

    count_samples(args.compiled_fasta_file)

if __name__ == "__main__":
    sys.exit(main())
