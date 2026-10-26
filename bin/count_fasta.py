#!/usr/bin/python3

import argparse
import logging
import sys
import csv

logging.basicConfig(level=logging.INFO, format="%(levelname)s : %(message)s")


class ReadReportParser(argparse.ArgumentParser):
    def error(self, msg):
        self.print_help()
        sys.stderr.write(f"\nERROR: {msg}\n")
        sys.exit(1)

if __name__ == "__main__":
    parser = ReadReportParser(
        description="Generate a CSV report with read counts"
    )

    parser.add_argument(
        "--sample_id",
        help="Sample ID for the report"
    )
    parser.add_argument(
        "--fasta",
        required=True,
        nargs="+",
        help="FASTA file"
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV file"
    )

    args = parser.parse_args()

FIELDS = [
    "id",
    "count"
    ]

def write_csv(path, row):
    """Write a single-row CSV."""
    try:
        with open(path, "w", newline="") as fpath:
            writer = csv.DictWriter(fpath, fieldnames=FIELDS)
            writer.writeheader()
            writer.writerow(row)
    except Exception as exp:
        logging.error(f"Cannot write CSV {path}: {exp}")
        sys.exit(1)

    logging.info(f"Report written: {path}")

def count_seqs(path):
    """Count sequences in FASTA."""

    try:
        logging.info(f"Counting sequences: {path}")
        seqs = len([1 for line in open(path) if line.startswith(">")])

    except Exception as e:
        logging.error(f"Cannot read FASTA {path}: {e}")
        sys.exit(1)

    return seqs

def process(fasta, output, sample_id):
    """Main processing function. """

    logging.info("Counting reads")
    fasta_counts = [count_seqs(f) for f in fasta]

    row = {
        "id": sample_id,
        "count": fasta_counts[0],
    }

    write_csv(output, row)

process(args.fasta, args.output, args.sample_id)