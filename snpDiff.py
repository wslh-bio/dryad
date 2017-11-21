#!/usr/bin/env python3
import subprocess as sub
import sys
import os
import argparse
from Bio import Align
from Bio import AlignIO
import csv

def getSNP(seq1,seq2):
    seq1 = seq1[0]
    seq2 = seq2[0]
    if len(seq1) == len(seq2):
        l = len(seq1)
        c = 0
        snp_count = 0
        while c < l:
            if (seq1[c] != seq2[c]) and (seq1[c] != '-') and (seq2[c] != '-'):
                #print("SNP found:"+str(seq1[c])+" -> "+str(seq2[c]))
                snp_count += 1
            c += 1
        return snp_count

#determine command line arguments and get path
parser = argparse.ArgumentParser(description='Annotate trees with SNP differences between branches')
parser.add_argument('-a',metavar='alignment', type=str,help="sequence alignment",required=True)
#parser.add_argument('-m',metavar='method', type=str,help="alignment method: clustal,emboss,fasta")
#parser.add_argument('-t',metavar='tree', type=str,help="newark tree",required=True)
args = parser.parse_args()
algn_path = os.path.abspath(args.a)
#tre_path = os.path.abspath(args.t)

#get alignment
align = AlignIO.read(algn_path,"fasta")

#initalize alignment and snp matrix
algn_matrix = []
snp_matrix = [[0]*len(align) for _ in range(len(align))]

#populate matrix dict
print("Populating data matrix")
for item in align:
    algn_matrix.append([item.id,list(item)])

#comparison
print("Running SNP comparison")
c = 0
while c < len(align):
    b = c + 1
    while b < len(align):
        print("comparing "+algn_matrix[c][0]+" to "+algn_matrix[b][0])
        snp_matrix[c][b] = getSNP(algn_matrix[c][1:],algn_matrix[b][1:])
        b += 1
    c += 1

with open('snp_table.csv','w') as csvfile:
    w = csv.writer(csvfile,delimiter=',',quotechar='|')
    w.writerow(['',list([item.id for item in align])])
    c = 0
    for row in zip(*snp_matrix):
        w.writerow([algn_matrix[c][0],list(row)])
        c += 1
