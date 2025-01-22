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
	parser.add_argument('version',
		help='Version of Dryad')

	return parser.parse_args(args)

def process_dfs(log, excluded, quast):

    df_log = pd.read_csv(log, sep='\t')
    df_excluded = pd.read_csv(excluded)

    if quast != 'empty.txt':
        df_quast = pd.read_csv(quast, sep='\t')
    else:
        df_quast = 'false'

    return df_log, df_excluded, df_quast

def join_dfs_no_quast(df_log, df_excluded, version):

    # Ensure sample is just name
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_excluded['Sample'] = df_excluded['Sample'].apply(lambda x: Path(x).stem)

    # Begin joining 
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')

    # Change float to string
    # Cluster Coverage (bps) excluded until this issue is resolved https://github.com/marbl/parsnp/issues/173
    df_log_excluded[['Sequence Length']] = df_log_excluded[['Sequence Length']].fillna(-1).astype(int)
   #  df_log_excluded[['Sequence Length','Cluster Coverage (bps)']] = df_log_excluded[['Sequence Length','Cluster Coverage (bps)']].fillna(-1).astype(int)

    # Change float to int and replace NA with -1
    df_log_excluded = df_log_excluded.replace(-1,'')

    # Rename columns
    df_log_excluded = df_log_excluded.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})

    # add Dryad version number
    df_log_excluded = df_log_excluded.assign(Version=version)

    df_log_excluded.to_csv('dryad_summary.csv', index=False)

def join_dfs_with_quast(df_log, df_excluded, df_quast, version):

    # Ensure sample is just name
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_excluded['Sample'] = df_excluded['Sample'].apply(lambda x: Path(x).stem)
    df_quast['Sample'] = df_quast['Sample'].apply(lambda x: Path(x).stem)

    # Begin joining 
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')
    df_log_excluded_quast = pd.merge(df_log_excluded, df_quast, on='Sample', how='outer')

    # Change float to int and replace NA with -1
    # Cluster Coverage (bps) excluded until this issue is resolved https://github.com/marbl/parsnp/issues/173
    df_log_excluded_quast[['Sequence Length','Contigs','N50']] = df_log_excluded_quast[['Sequence Length','Contigs','N50']].fillna(-1).astype(int)
    #df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']] = df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']].fillna(-1).astype(int)

    # Replace -1 with empty space
    df_log_excluded_quast = df_log_excluded_quast.replace(-1,'')

    # Rename and columns
    df_log_excluded_quast = df_log_excluded_quast.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})

    # add Dryad version number
    df_log_excluded_quast = df_log_excluded_quast.assign(Version=version)

    df_log_excluded_quast.to_csv('dryad_summary.csv', index=False)

def main(args=None):
    args = parse_args(args)

    l,e,q = process_dfs(args.aligner_log, args.excluded_samples, args.quast)

    if args.quast == 'empty.txt':
        join_dfs_no_quast(l,e,args.version)
    else:
        join_dfs_with_quast(l,e,q,args.version)

if __name__ == "__main__":
	sys.exit(main())