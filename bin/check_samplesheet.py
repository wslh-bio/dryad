#!/usr/bin/env python

import os
import sys
import errno
import argparse


def parse_args(args=None):
    Description = "Reformat input samplesheet file and check its contents."
    Epilog = "Example usage: python check_samplesheet.py <FILE_IN> <FILE_OUT>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FILE_IN", help="Input samplesheet file.")
    parser.add_argument("FILE_OUT", help="Output file.")

    return parser.parse_args(args)

# Creating a directory for the samples 
def make_dir(path):
    if len(path) > 0:
        try:
            os.makedirs(path)
        except OSError as exception:
            if exception.errno != errno.EEXIST:
                raise exception

# Print the error
def print_error(error, context="Line", context_str=""):
    error_str = "ERROR: Please check samplesheet -> {}".format(error)
    if context != "" and context_str != "":
        error_str = "ERROR: Please check samplesheet -> {}\n{}: '{}'".format(
            error, context.strip(), context_str.strip()
        )
    print(error_str)
    sys.exit(1)

# Checking the samplesheet to ensure it is formatted properly
def check_samplesheet(file_in, file_out):
    """
    This function checks that the samplesheet follows the following structure:

    sample,fasta
    SAMPLE_PE,SAMPLE_PE_RUN1.fa.gz
    SAMPLE_PE,SAMPLE_PE_RUN2.fa.gz
    SAMPLE_SE,SAMPLE_SE_RUN3.fa.gz
    """

    sample_mapping_dict = {}
    with open(file_in, "r") as fin:
        ## Check header
        MIN_COLS = 2
        HEADER = ["sample", "fasta"]
        header = [x.strip('"') for x in fin.readline().strip().split(",")]
        if header[: len(HEADER)] != HEADER:
            print("ERROR: Please check samplesheet header -> {} != {}".format(",".join(header), ",".join(HEADER)))
            sys.exit(1)

        ## Check sample entries
        for line in fin:
            cleaned_line = [x.strip().strip('"') for x in line.strip().split(",")]

            # Check valid number of columns per row
            if len(cleaned_line) < len(HEADER):
                print_error(
                    "Invalid number of columns (minimum = {})!".format(len(HEADER)),
                    "Line",
                    line,
                )
            num_cols = len([x for x in cleaned_line if x])
            if num_cols < MIN_COLS:
                print_error(
                    "Invalid number of populated columns (minimum = {})!".format(MIN_COLS),
                    "Line",
                    line,
                )

            ## Check sample name entries
            sample, fasta_1 = cleaned_line[: len(HEADER)] # Will look at # of elements in header
            if sample.find(" ") != -1:
                print(f"WARNING: Spaces have been replaced by underscores for sample: {sample}")
                sample = sample.replace(" ", "_")
            if not sample:
                print_error("Sample entry has not been specified!", "Line", line)

            ## Check FastA file extension
            for fasta in [fasta_1]:
                if fasta:
                    if fasta.find(" ") != -1:
                        print_error("FastA file contains spaces!", "Line", line)
                    if not fasta.endswith(".fasta.gz") and not fasta.endswith(".fa.gz") and not fasta.endswith(".fa"):
                        print_error(
                            "FastA file does not have extension '.fasta.gz', '.fa.gz', or '.fa'!",
                            "Line",
                            line,
                        )

            ## Auto-detect paired-end/single-end
            sample_info = []  ## [single_end, fasta_1]
            if sample and fasta_1 : 
                sample_info = [fasta_1]
            else:
                print_error("Invalid combination of columns provided!", "Line", line)

            ## Create sample mapping dictionary = { sample: [ fasta_1 ] }
            if sample not in sample_mapping_dict:
                sample_mapping_dict[sample] = [sample_info]
            else:
                if sample_info in sample_mapping_dict[sample]:
                    print_error("Samplesheet contains duplicate rows!", "Line", line)
                else:
                    sample_mapping_dict[sample].append(sample_info)

    ## Write validated samplesheet with appropriate columns
    if len(sample_mapping_dict) > 0:
        out_dir = os.path.dirname(file_out)

        make_dir(out_dir)

        with open(file_out, "w") as fout:
            fout.write(",".join(["sample", "fasta"]) + "\n")
            for sample in sorted(sample_mapping_dict.keys()):
                ## Check that multiple runs of the same sample are of the same datatype
                if not all(x[0] == sample_mapping_dict[sample][0][0] for x in sample_mapping_dict[sample]):
                    print_error(
                        "Multiple runs of a sample must be of the same datatype!",
                        "Sample: {}".format(sample),
                    )

                for idx, val in enumerate(sample_mapping_dict[sample]):
                    fout.write(",".join(["{}_T{}".format(sample, idx + 1)] + val) + "\n")
    else:
        print_error("No entries to process!", "Samplesheet: {}".format(file_in))

def main(args=None):
    args = parse_args(args)

    check_samplesheet(args.FILE_IN, args.FILE_OUT)
    print("done")

if __name__ == "__main__":
    sys.exit(main())