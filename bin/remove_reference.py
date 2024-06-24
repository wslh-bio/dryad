#!/usr/bin/env python3

import argparse
import logging
import os

from Bio.SeqRecord import SeqRecord
from Bio import SeqIO

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

parser =  argparse.ArgumentParser(description='Remove fasta from reference.')

parser.add_argument('compiled_fasta_file',
                    help='Complete fasta file that needs reference to be removed from.')
parser.add_argument('reference_name',
                    help='Original reference fasta name.')
args = parser.parse_args()

reference_name = os.path.basename(args.reference_name)
ref_fasta = "<" + reference_name + ".ref"

fasta_file = os.path.basename(args.compiled_fasta_file)
outFasta = "cleaned_" + fasta_file

with open(args.compiled_fasta_file, "r") as inFasta:
    print(args.compiled_fasta_file)
    for record in SeqIO.parse(inFasta, "fasta"):
      if record.id == ref_fasta:
         pass
      else: 
        SeqIO.write(record, outFasta, "fasta")

#remove fastafile.ref, any seq in fasta file seq record and then remove the record that has .ref in it