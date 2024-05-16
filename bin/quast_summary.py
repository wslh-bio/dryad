#!/usr/bin/python3.7
import os
import glob
import pandas as pd
from pandas import DataFrame

# function for summarizing quast output
def summarize_quast(file):
    # get sample id from file name and set up data list
    sample_id = os.path.basename(file).split(".")[0]
    # read in data frame from file
    df = pd.read_csv(file, sep='\\t')
    # get contigs, total length and assembly length columns
    df = df.iloc[:,[1,7,17]]
    # assign sample id as column
    df = df.assign(Sample=sample_id)
    # rename columns
    df = df.rename(columns={'# contigs (>= 0 bp)':'Contigs','Total length (>= 0 bp)':'Assembly Length (bp)'})
    # re-order data frame
    df = df[['Sample', 'Contigs','Assembly Length (bp)', 'N50']]
    return df

# get quast output files
files = glob.glob("data*/*.transposed.quast.report.tsv*")

# summarize quast output files
dfs = map(summarize_quast,files)
dfs = list(dfs)

# concatenate dfs and write data frame to file
if len(dfs) > 1:
    dfs_concat = pd.concat(dfs)
    dfs_concat.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
else:
    dfs = dfs[0]
    dfs.to_csv(f'quast_results.tsv',sep='\\t', index=False, header=True, na_rep='NaN')
