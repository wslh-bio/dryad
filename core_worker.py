#!/usr/bin/env python3
import subprocess as sub
import sys
import os
import docker
import time
from shutil import copyfile
#function to build core genome tree
def cg_tree(out,p_list,user_id,user_grp,client,threads):
    #keep track of progress
    #stage 1 - create temp dir
    #stage 2 - annotation
    #stage 3 - alignment
    #stage 4 - tree
    stage = 0

    if not os.path.isdir(out):
        sub.Popen(["mkdir",out]).wait()
    #create temp folder
    if os.path.isdir(out+'/dryad_temp'):
        stage = 1
    if stage == 0:
        sub.Popen(["mkdir",out+"/dryad_temp"]).wait()
        stage = 1
    oout = out
    out = out + '/dryad_temp'

    #examine temp file to get stage
    #temp file name
    temp_f = out+'/sc'
    try:
        with open(temp_f,'r') as st:
            check = st.read()
            if '2' in check:
                print('annotations were previously compleated skipping to alignment')
                stage = 2
            elif '3' in check:
                print('alignment previously compleated skipping to tree')
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
        print("annotating genomes in list")
        for gn in p_list:
            name = os.path.basename(gn)
            n = name.split('.')[0]
            copyfile(gn,out+'/'+name)
            client.containers.run("staphb/prokka","prokka --cpus {2} --outdir /data/{1} /data/{0}".format(name,n,threads),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
            print("compleated {}".format(name))

        #cleanup prokka output for next steps
        dirList = []
        ls = sub.Popen(["ls", out],stdout=sub.PIPE)
        dirList = ls.communicate()[0].decode('UTF-8').split()
        for item in dirList:
            if ".fasta" not in item:
                print('Copying '+item)
                sub.Popen("cp "+out+'/'+item+'/'+"*.gff"+' '+out+'/'+item+".gff",shell=True,stdout=sub.PIPE).wait()
        stage = 2
        with open(temp_f,'w') as st:
            st.write('2')

    if stage == 2:
        print("creating alignment")
        client.containers.run("staphb/roary","sh -c 'roary -e -mafft -p {0} -f roary_out *.gff'".format(threads),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)
        print("compleated alignment")

        #move core gene alignment out of roary_out
        sub.Popen(['cp',out+'/roary_out/core_gene_alignment.aln',out]).wait()
        stage = 3
        with open(temp_f,'w') as st:
            st.write('3')

    if stage == 3:
        print("creating maximum likelihood tree using 1000 bootstraps")
        client.containers.run("staphb/iqtree","iqtree -nt {0} -m GTR+G -bb 1000 -s core_gene_alignment.aln -pre cg_tree".format(threads),user=user_id+":"+user_grp, working_dir='/data', volumes={out:{'bind':'/data','mode':'rw'}}, remove=True)

        #naming based off time
        o_name = str(time.localtime().tm_year)[2:]+str(time.localtime().tm_mon)+str(time.localtime().tm_mday)+"_cg_tree.newick"

        if os.path.isfile(oout+'/'+o_name):
            c = 0
            while True:
                if not os.path.isfile(oout+'/'+o_name.split('.')[0] + "_{0}".format(c) + ".newick"):
                    o_name = o_name.split('.')[0] + "_{0}".format(c) + ".newick"
                    break
                c += 1

        #move tree out of temp folder
        print("writing out tree: {0}".format(o_name))
        sub.Popen(['cp',out+'/cg_tree.treefile',oout+'/'+o_name])
        stage = 4
        with open(temp_f,'w') as st:
            st.write('4')
