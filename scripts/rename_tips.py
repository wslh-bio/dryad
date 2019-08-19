#!/usr/bin/env python3
#rename_tips.py
#Kelsey Florek
#This script takes in a csv file and a tree file and uses the csv to find and replace the corresponding tip labels with the label from the csv [find,replace].
import subprocess as sub
import sys
import os
import argparse
import csv

#setup argparser to display help if no arguments
class MyParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write('error: %s\n' % message)
        self.print_help()
        sys.exit(2)

#determine command line arguments and get path
parser = MyParser(description='Find and replace tip labels using a csv: Column 1 - tip label, Column 2 - replacement tip label ')
parser.add_argument('csv', type=str,help="csv file used for find/replace (commma delimited)")
parser.add_argument('tree_file', type=str,help="Tree file to rename tips")

args = parser.parse_args()
csvfile = os.path.abspath(args.csv)
target = os.path.abspath(args.tree_file)

with open(csvfile,'r') as cfile:
    reader = csv.reader(cfile,delimiter=',')
    for row in reader:
        if row[0] != '' and row[1] != '':
            with open(target,'r+') as tfile:
                data = tfile.read()
                tfile.seek(0)
                output = ''
                for line in data.splitlines():
                    output += line.replace(row[0],row[1])+'\n'
                tfile.write(output)
                tfile.truncate()
