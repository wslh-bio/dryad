#! usr/bin/env python3

import pandas as pd
import argparse
import re
import sys

def parse_args(args=None):
	Description='Compares the samples input in the pipeline to the samples that exit parsnp'

	parser = argparse.ArgumentParser(description=Description)
	parser.add_argument('input_file',
		help='Samplesheet.valid.csv with all the samples denoted.')
	parser.add_argument('output_file',
		help='Parsnp tree output file to compare input samples to.')

	return parser.parse_args(args)

def create_list_of_input(samplesheet):

	samplesheet_df = pd.read_csv(samplesheet)
	all_samples = samplesheet_df['sample'].tolist()

	return all_samples

def check_parsnp_file(parsnp_file, input_list):

	with open(parsnp_file, "r") as infile:
		tree_content = infile.read()

	samples_from_tree = re.findall(r'[^:()]+', tree_content)

	# Clean up sample names, if necessary
	samples_from_tree = [sample.strip() for sample in samples_from_tree]

	# Step 3: Compare the lists
	missing_samples = set(input_list) - set(samples_from_tree)

	# Output results
	if missing_samples:
		print("Missing samples in parsnp.tree:")
		for sample in missing_samples:
			print(sample)
	else:
	    print("All samples are present in parsnp.tree.")

def main(args=None):
	args = parse_args(args)

	total_list = create_list_of_input(args.input_file)
	check_parsnp_file(args.output_file, total_list)

if __name__ == "__main__":
	sys.exit(main())