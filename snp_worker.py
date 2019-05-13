#!/usr/bin/env python3
import subprocess as sub
import sys
import os
import docker
import time
from shutil import copyfile
#function to build snp tree
def snp_tree(out,p_list,user_id,user_grp,client,reference,threads):
    lyveset_container = "staphb/lyveset:2.0.1"
    ninja_container = "nwflorek/ninja:1.2.2"
    #keep track of progress
    #stage 1 - create temp dir
    #stage 2 - trim and shuffel reads, set reference
    #stage 3 - build lyveset project
    #stage 4 - run lyveset
    stage = 0

    if not os.path.isdir(out):
        sub.Popen(["mkdir",out]).wait()
    #create temp folder
    if os.path.isdir(out+'/dryad_snp_temp'):
        stage = 1
    if stage == 0:
        sub.Popen(["mkdir",out+"/dryad_snp_temp"]).wait()
        stage = 1
    oout = out
    out = out + '/dryad_snp_temp'

    #start logging stdout and stderr
    #sys.stdout = open(out+'/'+'snp_'+str(os.getpid()) + ".out", "w",buffering=1)
    #sys.stderr = open(out+'/'+'snp_'+str(os.getpid()) + ".err", "w",buffering=1)

    #examine temp file to get stage
    #temp file name
    temp_f = out+'/ss'
    try:
        with open(temp_f,'r') as st:
            check = st.read()
            if '2' in check:
                print('reads were previously trimmed skipping to building project')
                stage = 2
            elif '3' in check:
                print('project previously built running Lyve-SET')
                reference = os.path.basename(reference)
                stage = 3
            elif '4' in check:
                print('alignment previously made, building tree')
                reference = os.path.basename(reference)
                stage = 4
            else:
                if stage == 1:
                    pass
                else:
                    print('unknown program status, delete temporary files and start again')
                    sys.exit()
    except FileNotFoundError:
        pass

    if stage == 1:
        print("trimming reads")
        for r1,r2 in zip(p_list[0::2],p_list[1::2]):
            #get id
            _id = os.path.basename(r1).split('_')[0]
            #rename raw reads
            r1raw = _id+"-raw_1.fastq.gz"
            r2raw = _id+"-raw_2.fastq.gz"
            #copy reads to temp dir
            copyfile(r1,out+'/'+r1raw)
            copyfile(r2,out+'/'+r2raw)
            #trim
            client.containers.run("staphb/trimmomatic","java -jar /Trimmomatic-0.38/trimmomatic-0.38.jar PE -threads {3} /data/{0} /data/{1} {2}_1.fastq.gz {2}_1U {2}_2.fastq.gz {2}_2U SLIDINGWINDOW:4:30".format(r1raw,r2raw,_id,threads),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
            #remove raw and unpaired reads
            sub.Popen(['rm',out+'/'+r1raw]).wait()
            sub.Popen(['rm',out+'/'+r2raw]).wait()
            sub.Popen(['rm',out+'/'+_id+'_1U']).wait()
            sub.Popen(['rm',out+'/'+_id+'_2U']).wait()
            #rename reads and compress for shuffeling
            print("compleated {0} / {1}".format(r1,r2))

        stage = 2
        with open(temp_f,'w') as st:
            st.write('2')

    if stage == 2:
        print("shuffle reads")
        client.containers.run(lyveset_container,"sh -c 'shuffleSplitReads.pl --numcpus {0} -o inter *.fastq.gz'".format(threads),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        print("compleated shuffle")
        print("building project")
        client.containers.run(lyveset_container,"set_manage.pl --create snp_tree",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        client.containers.run(lyveset_container,"sh -c 'for i in inter/*.fastq.gz; do set_manage.pl snp_tree --add-reads $i; done'",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        print("setting reference")
        sub.Popen(['cp',reference,out])
        reference = os.path.basename(reference)
        client.containers.run(lyveset_container,"sh -c 'set_manage.pl snp_tree --change-reference {0}'".format(reference),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        stage = 3
        with open(temp_f,'w') as st:
            st.write('3')

    if stage == 3:
        print("making SNP pseudo alignment with Lyve-SET")
        client.containers.run(lyveset_container,"sh -c 'launch_set.pl snp_tree --numcpus {0} --notrees -ref {1}'".format(threads,reference),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        stage = 4
        with open(temp_f,'w') as st:
            st.write('4')

    if stage == 4:
        #naming based off time
        out_prefix = str(time.localtime().tm_year)[2:]+str(time.localtime().tm_mon)+str(time.localtime().tm_mday) + ".nj"

        #incrament if already exists
        if os.path.isfile(oout+'/'+out_prefix+'.nj'):
            c = 0
            while True:
                if not os.path.isfile(oout+'/'+out_prefix.split('.')[0] + "_{0}".format(c) + ".nj"):
                    out_prefix = out_prefix.split('.')[0] + "_{0}".format(c) + ".nj"
                    break
                c += 1

        #build neghborhood joining tree
        client.containers.run(ninja_container,"sh -c 'phylip_converter.py snp_tree/msa/out.pairwiseMatrix.tsv phylip.pairwiseMatrix.tsv; ninja --in_type d phylip.pairwiseMatrix.tsv > {}'".format(out_prefix),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        #move matrix and alignment and tree out of temp folder
        sub.Popen(['cp',out+'/snp_tree/msa/out.informative.fasta',oout+'/'+out_prefix+'.aln.fna']).wait()
        sub.Popen(['cp',out+'/snp_tree/msa/out.pairwiseMatrix.tsv',oout+'/'+out_prefix+'.matrix.tsv']).wait()
        sub.Popen(['cp',out+'/'+out_prefix,oout+'/'+out_prefix]).wait()
