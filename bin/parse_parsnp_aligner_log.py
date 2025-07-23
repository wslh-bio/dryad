#!/usr/bin/python3.7
import re
import pandas as pd
import sys
import argparse
import logging

from pandas import DataFrame
from functools import reduce
from linecache import getline

# Reference for extracting lines after a match using enumerate and linecache:
# https://stackoverflow.com/questions/30286603/python-extract-text-4-lines-after-match

logging.basicConfig(level = logging.INFO, format = '%(levelname)s : %(message)s')

def parse_args(args=None):
    Description='Covert parsnpAligner.log to tsv'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('log',
                    help='Parsnp\'s parsnpAligner.log file.')
    parser.add_argument('add_reference',
                    help='Keep reference in output')
    return parser.parse_args(args)

def parse_log(log):
    logging.debug("Create empty dictionaries")
    lengthDict = {'Sample':[],'Sequence':[],'Sequence Length (bps)':[]}
    covDict = {'Sequence':[],'Cluster Coverage (%)':[]}

    with open(log,'r') as logFile:
        for ind, line in enumerate(logFile,1):
            logging.debug("Search for sequence lines")
            if re.search('Sequence \d+ : ',line):
                logging.debug("Get sequence # line and split at: 'Sequence # : Sample.fasta'")
                sampleLine = getline(logFile.name, ind).strip().split(":")
                logging.debug("Get length line (two lines after sequence line) and split at: 'Length: # bps'")
                lengthLine = getline(logFile.name, ind + 2).strip().split(":")
                logging.debug("Remove trailing white space from sequence #")
                sequence = sampleLine[0].rstrip()
                logging.debug("Remove leading white space from sample name")
                sample = sampleLine[1].lstrip()
                logging.debug("Remove leading white space and get rid of bps from length")
                length = lengthLine[1].lstrip().replace(" bps","")
                logging.debug("Append values to dictionary")
                lengthDict['Sample'].append(sample)
                lengthDict['Sequence'].append(sequence)
                lengthDict['Sequence Length (bps)'].append(length)
            logging.debug("Search for cluster coverage lines")
            if re.search('Cluster coverage in sequence \d+:',line):
                logging.debug("Split line at : 'Cluster coverage in sequence #:   #%'")
                sline = line.strip().split(":")
                logging.debug("Reformat line to Sequence # only")
                sequence = sline[0].replace('Cluster coverage in s','S')
                logging.debug("Remove % from cluster coverage percentage")
                clusterCov = sline[1].strip().replace('%','')
                logging.debug("Append values to dictionary")
                covDict['Sequence'].append(sequence)
                covDict['Cluster Coverage (%)'].append(clusterCov)
            logging.debug("Search for total coverage")
            if re.search('Total coverage among all sequences',line):
                totalCoverage = float(line.strip().split(":")[1].lstrip().replace('%',''))
    return lengthDict, covDict, totalCoverage

def createDF(lengthDict, covDict, totalCoverage, addRef):
    logging.debug("Convert length dictionary to data frame")
    lengthDF = pd.DataFrame.from_dict(lengthDict, orient='columns')
    logging.debug("Change sequence length from str to int")
    lengthDF['Sequence Length (bps)'] = lengthDF['Sequence Length (bps)'].astype(int)
    logging.debug("Convert coverage dictionary to data frame")
    covDF = pd.DataFrame.from_dict(covDict, orient='columns')
    logging.debug("Change cluster coverage from str to float")
    covDF['Cluster Coverage (%)'] = covDF['Cluster Coverage (%)'].astype(float)
    logging.debug("Add dfs to list")
    dfs = [lengthDF,covDF]
    logging.debug("Merge dfs in list")
    merged_df = reduce(lambda left,right: pd.merge(left,right,on=['Sequence'],how='left'), dfs)
    logging.debug("Calculate cluster coverage in bps")
    merged_df['Cluster Coverage (bps)'] = merged_df['Sequence Length (bps)'] * (merged_df['Cluster Coverage (%)']/100)
    merged_df['Cluster Coverage (bps)'] = merged_df['Cluster Coverage (bps)'].round().astype(int)
    logging.debug("Add total coverage column")
    merged_df = merged_df.assign(TotalCoverage=totalCoverage)
    merged_df.rename(columns={'TotalCoverage':'Total Coverage (%)'}, inplace = True)
    logging.debug("If reference is removed (default) drop reference from df")
    if addRef == "false":
        logging.debug("Drop reference row")
        merged_df.drop(merged_df[merged_df['Sample'].str.endswith('.ref')].index, inplace = True)
    logging.debug("Remove .contigs and .ref from sample names")
    merged_df['Sample'] = merged_df['Sample'].str.replace('.contigs', '')
    merged_df['Sample'] = merged_df['Sample'].str.replace('.ref', '')

    logging.debug("Write to file")
    merged_df.to_csv(f'aligner_log.tsv', sep='\t', index=False, header=True, na_rep='NaN')

def main(args=None):
    args = parse_args(args)
    lengthDict, covDict, totalCoverage = parse_log(args.log)
    createDF(lengthDict, covDict, totalCoverage, args.add_reference)

if __name__ == "__main__":
    sys.exit(main())
