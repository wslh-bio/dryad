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
            #copy reads to temp dir
            copyfile(r1,out+'/'+os.path.basename(r1))
            copyfile(r2,out+'/'+os.path.basename(r2))
            #move pointers to copied reads
            r1 = out+'/'+os.path.basename(r1)
            r2 = out+'/'+os.path.basename(r2)
            #unzip if compressed
            if '.gz' in r1:
                sub.Popen(['gunzip',r1])
            if '.gz' in r2:
                sub.Popen(['gunzip',r2])
            #set basename
            base_name = os.path.basename(r1).split('_')[0]
            #trim
            client.containers.run("nwflorek/trimassem","trimmomatic PE -threads 4 /data/{0} /data/{1} -baseout /data/{2} SLIDINGWINDOW:4:30".format(r1,r2,base_name),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
            #remove raw reads
            sub.Popen(['rm',r1])
            sub.Popen(['rm',r2])
            #rename reads for shuffeling
            sub.Popen(['mv',base_name+'_1P', base_name+'_1.fastq'])
            sub.Popen(['mv',base_name+'_2P', base_name+'_2.fastq'])
            print("compleated {0} / {1}".format(r1,r2))


        stage = 2
        with open(temp_f,'w') as st:
            st.write('2')

    if stage == 2:
        print("shuffle reads")
        client.containers.run("nwflorek/lyveset","shuffleSplitReads.pl --numcpus 4 -o inter *.fastq",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        print("compleated compleated shuffle")
        print("building project")
        client.containers.run("nwflorek/lyveset","set_manage.pl --create snp_tree",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        client.containers.run("nwflorek/lyveset","for i in inter/*.fastq.gz; do set_manage.pl snp_tree --add-reads $i done;",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        client.containers.run("nwflorek/lyveset","set_manage.pl snp_tree --change-reference {0}".format(reference),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        stage = 3
        with open(temp_f,'w') as st:
            st.write('3')

    if stage == 3:
        print("building SNP tree with Lyve-SET")
        client.containers.run("nwflorek/lyveset","launch_set.pl snp_tree --numcpus 4",user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

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
        sub.Popen(['cp',out+'snp_tree/msa/out.RAxML_bipartitions.raxml',oout+'/'+o_name])
        stage = 4
        with open(temp_f,'w') as st:
            st.write('4')

    if not keep_temp:
        print("remving temp dir: {0}".format(out))
        sub.Popen(["rm","-r",out]).wait()
