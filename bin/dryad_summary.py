#!/usr/bin/python3.7

import argparse
import sys
import logging

import pandas as pd

from pathlib import Path

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

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

    logging.debug("Creating dataframe from parsnp log and excluded log")
    df_log = pd.read_csv(log, sep='\t')
    df_excluded = pd.read_csv(excluded)

    logging.debug("Checks if there is or is not a quast file")
    if quast != 'empty.txt':
        df_quast = pd.read_csv(quast, sep='\t')
    else:
        df_quast = 'false'

    return df_log, df_excluded, df_quast

def join_dfs_no_quast(df_log, df_excluded, version):

    logging.debug("Setting the column order for no quast output")
    column_order = ["Sample","Excluded from Parsnp's analysis","Version"]

    logging.debug("Ensuring the sample is just the name")
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_excluded['Sample'] = df_excluded['Sample'].apply(lambda x: Path(x).stem)

    logging.debug("Merging the dataframes based on sample name")
    df_log_excluded = pd.merge(df_log, df_excluded, on='Sample', how='outer')

    logging.debug("Change float to string")
    logging.debug("Cluster Coverage (bps) excluded until this issue is resolved https://github.com/marbl/parsnp/issues/173")
    df_log_excluded[['Sequence Length']] = df_log_excluded[['Sequence Length']].fillna(-1).astype(int)
   #  df_log_excluded[['Sequence Length','Cluster Coverage (bps)']] = df_log_excluded[['Sequence Length','Cluster Coverage (bps)']].fillna(-1).astype(int)

    logging.debug("Change float to int and replace NA with -1")
    df_log_excluded = df_log_excluded.replace(-1,'')

    logging.debug("Rename columns")
    df_log_excluded = df_log_excluded.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})

    logging.debug("Add Dryad version number")
    df_log_excluded = df_log_excluded.assign(Version=version)

    logging.debug("Reordering columns based on column order")
    df_log_excluded = df_log_excluded.reindex(columns = column_order)

    logging.debug("Writing to csv")
    df_log_excluded.to_csv('dryad_summary.csv', index=False)

def join_dfs_with_quast(df_log, df_included, df_quast, version):

    logging.debug("Setting the column order for including quast output")
    column_order = ["Sample","Excluded from Parsnp's analysis","Contigs","N50","Assembly Length (bp)","Version"]

    logging.debug("Ensuring sample is just name")
    df_log['Sample'] = df_log['Sample'].apply(lambda x: Path(x).stem)
    df_included['Sample'] = df_included['Sample'].apply(lambda x: Path(x).stem)
    df_quast['Sample'] = df_quast['Sample'].apply(lambda x: Path(x).stem)

    logging.debug("Begin joining")
    df_log_included = pd.merge(df_log, df_included, on='Sample', how='outer')
    df_log_included_quast = pd.merge(df_log_included, df_quast, on='Sample', how='outer')

    logging.debug("Change float to int and replace NA with -1")
    logging.debug("Cluster Coverage (bps) excluded until this issue is resolved https://github.com/marbl/parsnp/issues/173")
    df_log_included_quast[['Sequence Length','Contigs','N50']] = df_log_included_quast[['Sequence Length','Contigs','N50']].fillna(-1).astype(int)
    #df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']] = df_log_excluded_quast[['Sequence Length','Cluster Coverage (bps)','Contigs','N50']].fillna(-1).astype(int)

    logging.debug("Replace -1 with empty space")
    df_log_included_quast = df_log_included_quast.replace(-1,'')

    logging.debug("Rename and columns")
    df_log_included_quast = df_log_included_quast.rename(columns={'excluded_from_analysis':'Excluded from Parsnp\'s analysis'})

    logging.debug("add Dryad version number")
    df_log_included_quast = df_log_included_quast.assign(Version=version)

    logging.debug("Reindexing based on column order")
    df_log_included_quast = df_log_included_quast.reindex(columns = column_order)

    logging.debug("Writing to csv")
    df_log_included_quast.to_csv('dryad_summary.csv', index=False)

def main(args=None):
    args = parse_args(args)

    l,e,q = process_dfs(args.aligner_log, args.excluded_samples, args.quast)

    if args.quast == 'empty.txt':
        join_dfs_no_quast(l,e,args.version)
    else:
        join_dfs_with_quast(l,e,q,args.version)

if __name__ == "__main__":
	sys.exit(main())