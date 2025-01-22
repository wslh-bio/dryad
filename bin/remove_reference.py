#!/usr/bin/env python3

import argparse
import os
import sys

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

def parse_args(args=None):
  Description='Remove fasta from reference.'

  parser = argparse.ArgumentParser(description=Description)
  parser.add_argument('compiled_fasta_file',
                      help='Complete fasta file that needs reference to be removed from.')
  return parser.parse_args(args)

def rename(compiled_fasta_file):
  fasta_file = os.path.basename(compiled_fasta_file)
  outFasta = "cleaned_" + fasta_file

  return outFasta

def process_file(compiled_fasta_file, outFasta):
  with open(compiled_fasta_file, "r") as inFasta, open(outFasta, "a") as output:
      for record in SeqIO.parse(inFasta, "fasta"):
        if record.id.endswith(".ref"):
           pass
        else: 
          SeqIO.write(record, output, "fasta")

def main(args=None):
    args = parse_args(args)
    outFasta = rename(args.compiled_fasta_file)
    process_file(args.compiled_fasta_file, outFasta)

if __name__ == "__main__":
    sys.exit(main())