#!/usr/bin/env python3

import sys
import os

def find_duplicates(f):
    isolates = []
    with open(f,'r') as read_names:
        for line in read_names:
            pnusa = os.path.basename(line)
            pnusa = pnusa.split('.')[0]
            if pnusa not in isolates:
                isolates.append(pnusa)
            else:
                print('found duplicate: \n',line,'\n')

path = sys.argv[1]
path = os.path.abspath(path)
find_duplicates(path)
