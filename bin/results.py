#!/usr/bin/python3.7

import os
import glob
import csv
import pandas as pd
from functools import reduce

with open('kraken_version.yml', 'r') as krakenFile:
    for l in krakenFile.readlines():
        if "kraken DB:" in l.strip():
            krakenDBVersion = l.strip().split(':')[1].strip()

files = glob.glob('*.tsv')
dfs = []

for file in files:
    df = pd.read_csv(file, header=0, delimiter='\\t')
    dfs.append(df)

# Merge tsvs frames
merged = reduce(lambda  left,right: pd.merge(left,right,on=['Sample'],how='left'), dfs)

# Merge comment columns and drop individual columns that were merged
cols = ['Comments', 'Comments_x', 'Comments_y']
merged['Combined'] = merged[cols].apply(lambda row: '; '.join(row.values.astype(str)), axis=1)
merged['Combined'] = merged['Combined'].str.replace('nan; ', '')
merged['Combined'] = merged['Combined'].str.replace('; nan', '')
merged['Combined'] = merged['Combined'].str.replace('contamination', 'Contamination')
merged['Combined'] = merged['Combined'].str.replace('nan', '')
merged.drop(cols,axis=1, inplace=True)

# Add kraken DB column
merged = merged.assign(krakenDB=krakenDBVersion)

# Add Workflow version column
merged = merged.assign(workflowVersion='${workflow.manifest.version}')

# Add NTC columns
ntc = merged[merged['Sample'].str.match('NTC')]

kraken_ntc_results = glob.glob("kraken_ntc_data/*")

ntc_result = 'PASS'
ntc_total_reads = []
ntc_SPN_reads = []
for file in kraken_ntc_results:
    id = file.split("/")[1].split(".kraken.txt")[0]
    spn_reads = 0
    total_reads = 0
    with open(file,'r') as csvfile:
        dialect = csv.Sniffer().sniff(csvfile.read(1024))
        csvfile.seek(0)
        reader = csv.reader(csvfile,dialect)
        for row in reader:
            if row[3] == "U":
                total_reads += int(row[1])
            if "root" in row[5]:
                total_reads += int(row[1])
            if row[4] == "1300":
                spn_reads += int(row[1])

    if total_reads >= ${params.ntc_read_limit}:
        ntc_result = "FAIL"
    if spn_reads >= ${params.ntc_spn_read_limit}:
        ntc_result = "FAIL"

ntc_total_reads.append(f"{id}: {total_reads}")
ntc_SPN_reads.append(f"{id}: {spn_reads}")

if not kraken_ntc_results:
    merged = merged.assign(ntc_reads="No NTC in data set")
    merged = merged.assign(ntc_spn="No NTC in data set")
    merged = merged.assign(ntc_result="FAIL")

else:
    merged = merged.assign(ntc_reads=", ".join(ntc_total_reads))
    merged = merged.assign(ntc_spn=", ".join(ntc_SPN_reads))
    merged = merged.assign(ntc_result=ntc_result)

merged = merged.rename(columns={'Contigs':'Contigs (#)','Combined':'Comments','ntc_reads':'Total NTC Reads','ntc_spn':'Total NTC SPN Reads','ntc_result':'NTC PASS/FAIL','krakenDB':'Kraken Database Version','workflowVersion':'SPNtypeID Version'})
merged = merged[['Sample','Contigs (#)','Assembly Length (bp)','N50','Median Coverage','Average Coverage','Pass Coverage','Total Reads','Reads Removed','Median Read Quality','Average Read Quality','Pass Average Read Quality','Percent Strep','Percent SPN', 'SecondGenus','Percent SecondGenus','Pass Kraken','Serotype','Comments','Kraken Database Version','SPNtypeID Version','Total NTC Reads','Total NTC SPN Reads','NTC PASS/FAIL']]
merged.to_csv('spntypeid_report.csv', index=False, sep=',', encoding='utf-8')
