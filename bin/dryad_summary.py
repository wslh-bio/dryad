#!/usr/bin/python3.7

import argparse
import sys

import pandas as pd

from pathlib import Path

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

def join_dfs_no_quast(df_log, df_excluded):

    # Ensure sample is just name
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_excluded['Sample'] = df_excluded['Sample'].apply(lambda x: Path(x).stem)

    # Begin joining 
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')

    # Change float to string
    df_log_excluded[['Sequence Length','Cluster Coverage (bps)']] = df_log_excluded[['Sequence Length','Cluster Coverage (bps)']].astype('string').str.split('.').str[0]

    # Rename and Drop columns
    df_log_excluded = df_log_excluded.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})

    df_log_excluded.to_csv('dryad_summary.csv', index=False)

def join_dfs_with_quast(df_log, df_excluded, df_quast):

    # Ensure sample is just name
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_excluded['Sample'] = df_excluded['Sample'].apply(lambda x: Path(x).stem)
    df_quast['Sample'] = df_quast['Sample'].apply(lambda x: Path(x).stem)

    # Begin joining 
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')
    df_log_excluded_quast = pd.merge(df_log_excluded, df_quast, on='Sample', how='outer')

    # Change float to int
    df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']] = df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']].fillna(-1).astype(int)

    # Rename and Drop columns
    df_log_excluded_quast = df_log_excluded_quast.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})
    df_log_excluded_quast = df_log_excluded_quast.drop('Assembly Length (bp)', axis=1)

    df_log_excluded_quast.to_csv('dryad_summary.csv', index=False)

def main(args=None):
    args = parse_args(args)

    l,e,q = process_dfs(args.aligner_log, args.excluded_samples, args.quast)

    if args.quast == 'false':
        join_dfs_no_quast(l,e)
    else:
        join_dfs_with_quast(l,e,q)

if __name__ == "__main__":
	sys.exit(main())