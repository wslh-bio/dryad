#!/usr/bin/python3.7

import pandas as pd
import argparse
import sys
import re
import logging

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
	Description='Compares the samples input in the pipeline to the samples that exit parsnp'
	Epilog='Usage: python3 compare_io.py <SAMPLESHEET> <LOG>'

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('input_file',
		help='Samplesheet.valid.csv with all the samples denoted.')
	parser.add_argument('output_file',
		help='Aligner log output file to compare input samples to.')

	return parser.parse_args(args)

def create_list_of_input(samplesheet):

	logging.debug("Setting up samplesheet dataframe which contains all samples.")
	samplesheet_df = pd.read_csv(samplesheet)
	all_samples = samplesheet_df['sample'].tolist()

	return all_samples

def check_parsnp_file(parsnp_log_file, input_list):

	logging.debug("Creating empty list to hold samples that were not present in parsnp file")
	not_present_list = []

	with open (parsnp_log_file, "r") as input:

		logging.debug("Reading contents of file")
		results = input.read()

		logging.debug("Going through each sample name in input list")
		for sample in input_list:

			logging.debug("Checking to see if sample is in contents of parsnp results")
			if sample in results:
				pass
			else:
				not_present_list.append(sample)

	return not_present_list

def main(args=None):
	args = parse_args(args)

	total_list = create_list_of_input(args.input_file)
	not_present_list = check_parsnp_file(args.output_file, total_list)

	df = pd.DataFrame({'Sample': total_list})
	df['excluded_from_analysis'] = df['Sample'].apply(lambda x: 'Yes' if x in not_present_list else 'No')
	df.to_csv('sample_exclusion_status.csv', index=False)

if __name__ == "__main__":
	sys.exit(main())