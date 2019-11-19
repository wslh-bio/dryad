#!/usr/bin/env python3
import sys, os

SRA_test_ids = ['SRR6484704','SRR7107020','SRR7155575','SRR7179892']

#grab SRA toolkit
print(f"\nDownloading SRA Toolkit")
os.system('wget http://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.10.0/sratoolkit.2.10.0-centos_linux64.tar.gz')
os.system('tar -xzf sratoolkit.2.10.0-centos_linux64.tar.gz')

#download dataset
print("\nDownloading the reads.")
for sra_id in SRA_test_ids:
    print(f"Downloading...{sra_id}")
    cmd = f'sratoolkit.2.10.0-centos_linux64/bin/fastq-dump --split-3 --gzip {sra_id}'
    os.system(cmd)

#clean up
os.system('rm -rf sratoolkit.2.10.0-centos_linux64')
os.remove('sratoolkit.2.10.0-centos_linux64.tar.gz')

#fetch reference sequence
print(f"\nDownloading reference sequence.")
os.system('wget http://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz')
os.system('gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz')

#create full read list
with open('test_read_list.txt','w') as outfile:
    for sra_id in SRA_test_ids:
        outfile.write(f'{sra_id}_1.fastq.gz\n')
        outfile.write(f'{sra_id}_2.fastq.gz\n')

#run both core genome and snp
os.system('./dryad')
os.system('./dryad cg')
os.system('./dryad snp')
os.system('./dryad all test_read_list.txt GCF_000005845.2_ASM584v2_genomic.fna -t 2')
