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
	parser.add_argument('--quast_results',
		help='Supplies quast file, if run')
	parser.add_argument('--run_name',
        type=str,
		help='This is supplied by the nextflow config and can be changed via the usual methods i.e. command line.'),
	parser.add_argument('--dryad_version',
		help='Version of Dryad')

	return parser.parse_args(args)

def create_report(df_quast, version, WFRunName):

    logging.debug("Setting the column order for including quast output")
    column_order = ['Sample',
                    'Contigs',
                    'N50',
                    'Assembly Length (bp)',
                    'Version']

    logging.debug("Ensuring sample is just name")
    df_quast['Sample'] = df_quast['Sample'].apply(lambda x: Path(x).stem)

    logging.debug("Change float to int and replace NA with -1")
    df_quast[['Contigs',
              'N50']] = df_quast[['Contigs',
                                  'N50']].fillna(-1).astype(int)

    logging.debug("Replace -1 with empty space")
    df_quast = df_quast.replace(-1,'')

    logging.debug("Add Dryad version number")
    df_quast = df_quast.assign(Version=version)

    logging.debug("Reorder columns based on column order")
    df_quast = df_quast.reindex(columns = column_order)

    logging.debug("Writing to csv")
    df_quast.to_csv(f'{WFRunName}_dryad_summary.csv', index=False)

def main(args=None):
    args = parse_args(args)

    logging.debug("Creating dataframe from quast results")
    q = pd.read_csv(args.quast_results, sep='\t')

    create_report(q,args.dryad_version,args.run_name)

if __name__ == "__main__":
	sys.exit(main())