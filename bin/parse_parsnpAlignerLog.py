import re
from linecache import getline
import pandas as pd
from pandas import DataFrame
from functools import reduce
import sys
import argparse

# Reference for extracting lines after a match using enumerate and linecache:
# https://stackoverflow.com/questions/30286603/python-extract-text-4-lines-after-match

def parse_args(args=None):
    Description='Covert parsnpAligner.log to tsv'

    parser = argparse.ArgumentParser(description=Description)
    parser.add_argument('log',
                    help='Parsnp\'s parsnpAligner.log file.')
    parser.add_argument('add_reference',
                    help='Keep reference in output')
    return parser.parse_args(args)

def parse_log(log):
    # create empty dictionaries
    lengthDict = {'Sample':[],'Sequence':[],'Sequence Length':[]}
    covDict = {'Sequence':[],'Cluster Coverage (%)':[]}

    with open(log,'r') as logFile:
        for ind, line in enumerate(logFile,1):
            # search for sequence lines
            if re.search('Sequence \d+ : ',line):
                # get sequence # line and split at : (sequence on the left and sample on the right)
                # Sequence # : Sample.fasta
                sampleLine = getline(logFile.name, ind).strip().split(":")
                # get length line (two lines after sequence line) and split at : (sequence length on the right)
                # Length:   # bps
                lengthLine = getline(logFile.name, ind + 2).strip().split(":")
                # remove trailing white space from sequence #
                sequence = sampleLine[0].rstrip()
                # remove leading white space from sample name
                sample = sampleLine[1].lstrip()
                # remove leading white space and get rid of bps from length
                length = lengthLine[1].lstrip().replace(" bps","")
                # append values to dictionary
                lengthDict['Sample'].append(sample)
                lengthDict['Sequence'].append(sequence)
                lengthDict['Sequence Length'].append(length)
            # search for cluster coverage lines
            if re.search('Cluster coverage in sequence \d+:',line):
                # split line at :
                # Cluster coverage in sequence #:   #%
                sline = line.strip().split(":")
                # reformat line to Sequence # only
                sequence = sline[0].replace('Cluster coverage in s','S')
                # remove % from cluster coverage percentage
                clusterCov = sline[1].strip().replace('%','')
                # append values to dictionary
                covDict['Sequence'].append(sequence)
                covDict['Cluster Coverage (%)'].append(clusterCov)
            # search for total coverage
            if re.search('Total coverage among all sequences',line):
                totalCoverage = float(line.strip().split(":")[1].lstrip().replace('%',''))
    return lengthDict, covDict, totalCoverage

def createDF(lengthDict, covDict, totalCoverage, addRef):
    # convert length dictionary to data frame
    lengthDF = pd.DataFrame.from_dict(lengthDict, orient='columns')
    # change sequence length from str to int
    lengthDF['Sequence Length'] = lengthDF['Sequence Length'].astype(int)

    # convert coverage dictionary to data frame
    covDF = pd.DataFrame.from_dict(covDict, orient='columns')
    # change cluster coverage from str to float
    covDF['Cluster Coverage (%)'] = covDF['Cluster Coverage (%)'].astype(float)

    # add dfs to list
    dfs = [lengthDF,covDF]

    # merge dfs in list
    merged_df = reduce(lambda left,right: pd.merge(left,right,on=['Sequence'],how='left'), dfs)

    # calculate cluster coverage in bps
    merged_df['Cluster Coverage (bps)'] = merged_df['Sequence Length'] * (merged_df['Cluster Coverage (%)']/100)
    merged_df['Cluster Coverage (bps)'] = merged_df['Cluster Coverage (bps)'].round().astype(int)

    # add total coverage column
    merged_df = merged_df.assign(TotalCoverage=totalCoverage)
    merged_df.rename(columns={'TotalCoverage':'Total Coverage (%)'}, inplace = True)

    # if reference is removed (default) drop reference from df 
    if addRef == FALSE:
        # drop reference row
        merged_df.drop(merged_df[merged_df['Sample'].str.endswith('.ref')].index, inplace = True)

    # write to file
    merged_df.to_csv(f'aligner_log.tsv', sep='\t', index=False, header=True, na_rep='NaN')

def main(args=None):
    args = parse_args(args)
    lengthDict, covDict, totalCoverage = parse_log(args.log)
    createDF(lengthDict, covDict, totalCoverage, args.add_reference)

if __name__ == "__main__":
    sys.exit(main())
