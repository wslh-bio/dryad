#!/bin/bash
#This script impliments use of the lyve-SET SNP tree construction pipline on the CHTC condor grid.
#
#Order of arguments:
#1. path to reads
#2. path to reference
#3. number of cpus
#4. project name

# Shuffle your reads if they are not already. This command
# creates a folder interleaved and creates interleaved files
shuffleSplitReads.pl --numcpus $3 -o interleaved $1/*.fastq.gz
# Create the project directory
set_manage.pl --create $4
# Add reads
for i in interleaved/*.fastq.gz; do
  set_manage.pl $4 --add-reads $i
done;
# Specify your reference genome
set_manage.pl $4 --change-reference $2

#launch lyve-SET
launch_set.pl $4 --numcpus $3
