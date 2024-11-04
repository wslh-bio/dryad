#!/usr/bin/python3.7

import argparse
import sys

import pandas as pd

def parse_args(args=None):
	Description='Summarized both alignment free and alignment based output from Dryad.'
	Epilog='Usage: python3 dryad_summary.py '

	parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
	parser.add_argument('aligner_log',
		help='Output of parse_parsnp_aligner_log.py')
	parser.add_argument('quast',
		help='Supplies quast file, if run')
	parser.add_argument('excluded_samples',
		help='Output of compare_io.py')

	return parser.parse_args(args)

def process_dfs(log, excluded, quast):

    df_log = pd.read_csv(log, sep='\t')
    df_excluded = pd.read_csv(excluded)

    if quast != 'false':
        df_quast = pd.read_csv(quast, sep='\t')
    else:
        df_quast = 'false'

    return df_log, df_excluded, df_quast

def join_dfs(df_log, df_excluded, df_quast):
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')
    print(df_log_excluded)

def main(args=None):
    args = parse_args(args)

    process_dfs(args.aligner_log, args.excluded_samples, args.quast)

if __name__ == "__main__":
	sys.exit(main())