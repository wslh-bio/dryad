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
  parser.add_argument('reference_name',
                      help='Original reference fasta name.')
  return parser.parse_args(args)

def rename(reference_name, compiled_fasta_file):
  reference_name = os.path.basename(reference_name)
  ref_fasta = reference_name + ".ref"

  fasta_file = os.path.basename(compiled_fasta_file)
  outFasta = "cleaned_" + fasta_file

  return ref_fasta, outFasta

def process_file(compiled_fasta_file, ref_fasta, outFasta):
  with open(compiled_fasta_file, "r") as inFasta, open(outFasta, "a") as output:
      for record in SeqIO.parse(inFasta, "fasta"):
        if record.id == ref_fasta:
           pass
        else: 
          SeqIO.write(record, output, "fasta")

def main(args=None):
    args = parse_args(args)

    ref_fasta, outFasta = rename(args.reference_name, args.compiled_fasta_file)
    process_file(args.compiled_fasta_file, ref_fasta, outFasta)

if __name__ == "__main__":
    sys.exit(main())
