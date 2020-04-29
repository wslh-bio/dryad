#!/usr/bin/env python3
#Author: Kelsey Florek and Abigail Shockey
#email: kelsey.florek@slh.wisc.edu, abigail.shockey@slh.wisc.edu
#description: A pipeline for constructing SNP based and  core gene set reference free phylogenies


import sys,os,re
import argparse
from shutil import which, copyfile
from datetime import date
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

    parser = MyParser(description='Rebuild a previously generated PDF report.')
    parser.add_argument('--snp_matrix',type=str,help="path to snp matrix")
    parser.add_argument('--cg_tree',type=str,help="path to core genome tree")
    parser.add_argument('--ar',type=str,help="path to ar TSV file")
    parser.add_argument('--rmd',type=str,help="path to Rmarkdown file (.Rmd)")
    parser.add_argument('--config','-c', type=str,help="Nextflow custom configureation")
    parser.add_argument('--get_config',action="store_true",help="get a Nextflow configuration template for dryad")

    args = parser.parse_args()

    #give config to user if requested
    if args.get_config:
        config_path = os.path.join(dryad_path,"configs/dryad_config_template.config")
        dest_path = os.path.join(os.getcwd(),date.today().strftime("%y-%m-%d")+"_dryad.config")
        copyfile(config_path,dest_path)
        sys.exit()

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

    if args.rebuild_report:
        snp_mat = os.path.abspath(args.snp_matrix)
        cg_tree = os.path.abspath(args.cg_tree)
        ar_tsv = os.path.abspath(args.ar)
        rmd = os.path.abspath(args.rmd)

        #build command
        command = nextflow_path
        command = command + f" {config} run {dryad_path}/rebuild_report.nf {profile} --snp_matrix {snp_mat} --cg_tree {cg_tree} --ar_tsv {ar_tsv} --rmd {rmd} {work}"

    #run command using nextflow in a subprocess
    print("Starting the Dryad pipeline:")
    child = pexpect.spawn(command)
    child.interact()
