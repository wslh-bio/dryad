#! usr/bin/env python3

import pandas as pd
import argparse
import sys

def parse_args(args=None):
	Description='Compares the amount of samples input to the pipeline to the samples that exit parsnp'

	parser = argparse.ArgumentParser(description=Description)
	parser.add_argument('input_file',
		help='Samplesheet.valid.csv with all the samples denoted.')
	parser.add_argument('output_file',
		help='Parsnp output file to compare input samples to.')

	return parser.parse_args(args)

def create_list_of_input(samplesheet):
	input_list = []
	count = 0

	with open(samplesheet, "r") as infile:
		for line in infile:
			input_list.append(line.split(",")[0])
			count +=1

	return input_list, count

def check_parsnp_file(parsnp_file, input_list):
	parsnp_list = []
	missing_list = []
	with open(parsnp_file, "r") as infile:
		for line in infile:
			print(line + '\n')
			x = line.split(":")
			print(x)
			
			if x in input_list:
				parsnp_list.append(x)
			else:
				if type(x) == 'int':
					pass
				else:
					missing_list.append(x)

def main(args=None):
	args = parse_args(args)

	total_list, total_count = create_list_of_input(args.input_file)
	check_parsnp_file(args.output_file, total_list)

if __name__ == "__main__":
	sys.exit(main())