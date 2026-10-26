#!/usr/bin/env python3

import os
import sys
import glob
import argparse
import logging

import pandas as pd
from functools import partial

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def summarize_quast(file):
    logging.debug("Get sample id from file name and set up data list")
    sample_id = os.path.basename(file).split('.')[0]

    logging.debug("Read in data frame from file")
    df = pd.read_csv(file, sep='\t')

    logging.debug("Get contigs, total length and assembly length columns")
    df = df.loc[:,['# contigs','Total length', 'N50']]

    logging.debug("Assign sample id as column")
    df = df.assign(Sample=sample_id)

    logging.debug("Rename columns")
    df = df.rename(columns={'# contigs':'Contigs','Total length':'Assembly Length (bp)'})

    logging.debug("Re-order data frame")
    df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]

    return df

def grab_files():

    logging.info("Obtaining all QUAST output files")
    files = glob.glob('data*/*.transposed.quast.report.tsv*')

    return files

def summarize_output(files):

    summarize_quast_partial = partial(summarize_quast)

    logging.info("Summarizing quast output files")
    dfs = map(summarize_quast_partial,files)
    dfs = list(dfs)

    return dfs

def concatenate_dfs(dfs):
    logging.debug("Concatenate dfs and write data frame to file")
    if len(dfs) > 1:
        dfs_concat = pd.concat(dfs)
        dfs_concat.to_csv('quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')
    else:
        dfs = dfs[0]
        dfs.to_csv('quast_results.tsv',sep='\t', index=False, header=True, na_rep='NaN')

logging.info("Begin compiling all results for output file.")
files = grab_files()
dfs = summarize_output(files)
concatenate_dfs(dfs)