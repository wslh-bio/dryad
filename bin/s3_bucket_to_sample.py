#!/usr/bin/env python

import os
import sys
import argparse
import boto3

def parse_args(args=None):
    Description = (
        "Generate Dryad samplesheet from a S3 bucket of FastA files."
    )
    Epilog = "Example usage: python s3_bucket_to_samplesheet.py <FASTA_S3> <SAMPLESHEET_FILE>"

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument("FASTA_S3", help="S3 location (S3://bucket-name/key-name) containing raw FastA files.")
    parser.add_argument("SAMPLESHEET_FILE", help="Output samplesheet file.")
    parser.add_argument(
        "-a",
        "--assembly_extension",
        type=str,
        dest="ASSEMBLY_EXTENSION",
        default=".fa",
        help="File extension for assembly.",
    )
    parser.add_argument(
        "-sn",
        "--sanitise_name",
        dest="SANITISE_NAME",
        action="store_true",
        help="Whether to further sanitise FastA file name to get sample id. Used in conjunction with --sanitise_name_delimiter and --sanitise_name_index.",
    )
    parser.add_argument(
        "-sd",
        "--sanitise_name_delimiter",
        type=str,
        dest="SANITISE_NAME_DELIMITER",
        default="_",
        help="Delimiter to use to sanitise sample name.",
    )
    parser.add_argument(
        "-si",
        "--sanitise_name_index",
        type=int,
        dest="SANITISE_NAME_INDEX",
        default=1,
        help="After splitting FastA file name by --sanitise_name_delimiter all elements before this index (1-based) will be joined to create final sample name.",
    )
    return parser.parse_args(args)


def fasta_s3_to_samplesheet(
    fasta_s3,
    samplesheet_file,
    assembly_extension=".fa",
    sanitise_name=False,
    sanitise_name_delimiter="_",
    sanitise_name_index=1,
):
    ## Initialize AWS Client
    s3_client = boto3.client('s3')

    ##s Get bucket and prefix
    bucket_prefix = fasta_s3.strip('s3://').strip('S3://').split('/')
    bucket = bucket_prefix[0]
    prefix = '/'.join(bucket_prefix[1:])

    def sanitize_sample(path, extension):
        """Retrieve sample id from filename"""
        sample = os.path.basename(path).replace(extension, "")
        if sanitise_name:
            sample = sanitise_name_delimiter.join(
                os.path.basename(path).split(sanitise_name_delimiter)[
                    :sanitise_name_index
                ]
            )
        return sample

    def get_fastas(extension):
        ## object holder
        object_list = []

        ## get objects from s3 and append to list
        response = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix)
        for content in response['Contents']:
            if extension in content['Key']:
                object_list.append(os.path.join('s3://',response['Name'],content['Key']))
    
        ## if content is in parts fetch other parts
        if response['IsTruncated']:
            while(response['IsTruncated']):
                response = s3_client.list_objects_v2(Bucket=bucket, Prefix=prefix, ContinuationToken=response['NextContinuationToken'])
                for content in response['Contents']:
                    if extension in content['Key']:
                        object_list.append(os.path.join('s3://',response['Name'],content['Key']))

        return sorted(object_list)

    read_dict = {}

    ## Get assembly files
    for assembly_file in get_fastas(assembly_extension):
        sample = sanitize_sample(assembly_file, assembly_extension)
        if sample not in read_dict:
            read_dict[sample] = {"assembly": []}
        read_dict[sample]["assembly"].append(assembly_file)

    ## Write to file
    if len(read_dict) > 0:
        out_dir = os.path.dirname(samplesheet_file)
        if out_dir and not os.path.exists(out_dir):
            os.makedirs(out_dir)

        with open(samplesheet_file, "w") as fout:
            header = ["sample", "fasta"]
            fout.write(",".join(header) + "\n")
            for sample, reads in sorted(read_dict.items()):
                for idx, assembly in enumerate(reads["assembly"]):
                    sample_info = ",".join([sample, assembly])
                    fout.write(f"{sample_info}\n")
    else:
        error_str = (
            "\nWARNING: No FastA files found so samplesheet has not been created!\n\n"
        )
        error_str += "Please check the values provided for the:\n"
        error_str += "  - Path to the directory containing the FastA files\n"
        error_str += "  - '--assembly_extension' parameter\n"
        print(error_str)
        sys.exit(1)


def main(args=None):
    args = parse_args(args)

    fasta_s3_to_samplesheet(
        fasta_s3=args.FASTA_S3,
        samplesheet_file=args.SAMPLESHEET_FILE,
        assembly_extension=args.ASSEMBLY_EXTENSION,
        sanitise_name=args.SANITISE_NAME,
        sanitise_name_delimiter=args.SANITISE_NAME_DELIMITER,
        sanitise_name_index=args.SANITISE_NAME_INDEX,
    )


if __name__ == "__main__":
    sys.exit(main())