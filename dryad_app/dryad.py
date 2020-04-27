#!/usr/bin/env python3
#Author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#description: A pipeline for constructing SNP based and  core gene set reference free phylogenies


import sys,os,re
import argparse
from shutil import which
import pexpect

import re, sys

def main():
    #get nextflow executable
    lib_path = os.path.abspath(os.path.dirname(__file__) + '/' + '../lib')
    dryad_path = os.path.abspath(os.path.dirname(__file__))
    nextflow_path = os.path.join(lib_path,'nextflow')

    #setup argparser to display help if no arguments
    class MyParser(argparse.ArgumentParser):
        def error(self, message):
            sys.stderr.write('error: %s\n' % message)
            self.print_help()
            sys.exit(2)

    parser = MyParser(description='A comprehensive tree building program.')
    parser.add_argument('reads_path', type=str,help="Path to the location of the raw reads in the fastq format.")
    parser.add_argument('--output','-o',metavar="<output_path>",type=str,help="Path to ouput directory, default \"dryad_results\".",default="dryad_results")
    parser.add_argument('--core-genome','-cg',default=False, action="store_true", help="Construct a core-genome tree.")
    parser.add_argument('--snp','-s',default=False, action="store_true", help="Construct a SNP tree. Note: Requires a reference genome in fasta format (-r).")
    parser.add_argument('-r',metavar='<path>', type=str,help="Reference genome for SNP pipeline.")
    parser.add_argument('-ar',default=False, action="store_true", help="Detect AR mechanisms.")
    parser.add_argument('--sep',metavar="sep_chars",type=str,help="Dryad identifies sample names from the name of the read file by splitting the name on the specified separating characters, default \"_\".",default="_")
    parser.add_argument('--profile', type=str,choices=["docker", "singularity"],help="Nextflow profile. Default will try docker first, then singularity if the docker executable cannot be found.")
    parser.add_argument('--config','-c', type=str,help="Nextflow custom configureation.")
    parser.add_argument('--get_config',action="store_true",help="Get a Nextflow configuration template for dryad.")
    parser.add_argument('--resume', default="", action="store_const",const="-resume",help="resume a previous run")


    args = parser.parse_args()

    #give config to user if requested
    if args.get_config:
        config_path = os.path.join(dryad_path,"configs/dryad_config_template.config")
        dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_dryad.config")
        copyfile(config_path,dest_path)
        sys.exit()

    #check for reference sequence
    if args.snp and args.r == None:
        parser_dryad.print_help()
        print("Please specify a reference sequence for the SNP pipeline.")
        sys.exit(1)

    #check if we are using docker or singularity
    if which('docker'):
        profile = '-profile docker'
    elif which('singularity'):
        profile = '-profile singularity'
    else:
        profile = ''

    #check for config or profile
    config = ""
    if args.config:
        config = "-C " + os.path.abspath(args.config)
        profile = ""
    elif args.profile:
        profile = args.profile
    elif not profile:
        print('Singularity or Docker is not installed or not in found in PATH.')
        sys.exit(1)

    #set work dir into local logs dir if profile not aws
    work = ""
    if profile:
        work = f"-w {args.output}/logs/work"

    #build nextflow command
    selections = ""
    if args.ar:
        selections += " --ar"
    if args.core_genome:
        selections += " --cg"
    if args.snp:
        selections += f" --snp --snp_reference {args.r}"
    #add other arguments
    other_args = f"--name_split_on {args.sep} --outdir {args.output}"
    #build command
    command = nextflow_path
    command = command + f" {config} run {dryad_path}/dryad.nf {profile} {args.resume} --reads {args.reads_path} {selections} {other_args} {work}"

    #run command using nextflow in a subprocess
    print("Starting the Dryad pipeline:")
    child = pexpect.spawn(command)
    child.interact()
