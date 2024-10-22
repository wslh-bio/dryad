#!/usr/bin/python3.7

import pandas as pd
import argparse
import sys

from pathlib import Path

def parse_args(args=None):
	Description='Compares the samples input in the pipeline to the samples that exit parsnp'
	Epilog='Usage: python3 compare_io.py <SAMPLESHEET> <PARSNP.TREE>'

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('input_file',
		help='Samplesheet.valid.csv with all the samples denoted.')
	parser.add_argument('output_file',
		help='Aligner log output file to compare input samples to.')

	return parser.parse_args(args)

def create_list_of_input(samplesheet):

	samplesheet_df = pd.read_csv(samplesheet)
	all_samples = samplesheet_df['sample'].tolist()

	return all_samples

def check_parsnp_file(parsnp_file, input_list):

	# Setting up empty lists
	not_present_list = []
	cleaned_samples = []

	parsnp_df = pd.read_csv(parsnp_file,sep='\t')
	samples_present = parsnp_df.loc[:'Sample']

	# Clean up sample names, if necessary
	for sample in samples_present:
		sample = Path(sample).stem 
		cleaned_samples.append(sample)

	# Compare the lists
	for sample in input_list:
		if sample not in cleaned_samples:
			not_present_list.append(sample)

	return not_present_list

def main(args=None):
	args = parse_args(args)

	total_list = create_list_of_input(args.input_file)
	not_present_list = check_parsnp_file(args.output_file, total_list)

	with open("excluded_samples_from_parsnp.txt", 'w') as infile:

		if len(not_present_list) == 0:
			infile.write("No samples were excluded from parsnp's analysis.")

		else:
			infile.write("Samples were excluded from parsnp's analysis.\n")
			for sample in not_present_list:
				infile.write(f"{sample}\n")

if __name__ == "__main__":
	sys.exit(main())