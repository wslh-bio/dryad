#!/usr/bin/env python3
import subprocess as sub
import sys
import os
import docker
import time
from shutil import copyfile
#function to build snp tree
def snp_tree(out,p_list,user_id,user_grp,client,reference,keep_temp):
    #keep track of progress
    #stage 1 - create temp dir
    #stage 2 - trim reads
    #stage 3 - shuffle and build lyveset project
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
                stage = 3
            elif '4' in check:
                print('tree previously completed cleaning up')
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
            client.containers.run("nwflorek/trimassem","trimmomatic PE -threads 4 /data/{0} /data/{1} {2}_1.fastq.gz {2}_1U {2}_2.fastq.gz {2}_2U SLIDINGWINDOW:4:30".format(r1raw,r2raw,_id),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
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
        client.containers.run("nwflorek/lyveset","sh -c 'shuffleSplitReads.pl --numcpus 4 -o inter *.fastq.gz'",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        print("compleated shuffle")
        print("building project")
        client.containers.run("nwflorek/lyveset","set_manage.pl --create snp_tree",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        client.containers.run("nwflorek/lyveset","sh -c 'for i in inter/*.fastq.gz; do set_manage.pl snp_tree --add-reads $i; done'",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        stage = 3
        with open(temp_f,'w') as st:
            st.write('3')

    if stage == 3:
        sub.Popen(['cp',reference,out])
        reference = os.path.basename(reference)
        print("building SNP tree with Lyve-SET")
        client.containers.run("nwflorek/lyveset","sh -c 'set_manage.pl snp_tree --change-reference {0}'".format(reference),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        client.containers.run("nwflorek/lyveset","sh -c 'launch_set.pl snp_tree --numcpus 1'",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        #naming based off time
        o_name = str(time.localtime().tm_year)[2:]+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+"_snp_tree.raxml"

        if os.path.isfile(oout+'/'+o_name):
            c = 0
            while True:
                if not os.path.isfile(oout+'/'+o_name.split('.')[0] + "_{0}".format(c) + ".raxml"):
                    o_name = o_name.split('.')[0] + "_{0}".format(c) + ".raxml"
                    break
                c += 1

        #move tree out of temp folder
        print("writing out tree: {0}".format(o_name))
        sub.Popen(['cp',out+'/snp_tree/msa/out.RAxML_bipartitions.raxml',oout+'/'+o_name]).wait()
        stage = 4
        with open(temp_f,'w') as st:
            st.write('4')

    if not keep_temp:
        print("remving temp dir: {0}".format(out))
        sub.Popen(["rm","-r",out]).wait()
