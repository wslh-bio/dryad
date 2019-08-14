#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#core genome pipeline


import os,sys
import multiprocessing as mp
import psutil
import shutil

from app.lib import getfiles,checkexists,check_update_status
import app.calldocker as cd


#assembly function
def assemble_reads(jobs,cpu_job,outdir):
    logfile = os.path.join(outdir,'assembly.log')
    input_path = os.path.join(outdir,'trimmed')

    #determine free ram
    free_ram = int(psutil.virtual_memory()[1]/1000000000)
    ram_job = int(free_ram / jobs)

    #assemble
    assemblies_path = os.path.join(outdir,"assemblies")
    checkexists(assemblies_path)

    #get trimmed reads
    fastqs = list(getfiles(input_path))[0]
    cmds = []
    read_path = ''
    for read_pair in fastqs:
        #main command
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])
        sid = os.path.basename(read_pair[0]).split('_')[0]

        cmds.append('bash -c \"shovill --R1 /data/{0} --R2 /data/{1} --outdir /output/{2} --cpus {3} --ram {4} && mv /output/{2}/contigs.fa /output/{2}/{2}.fa \"'.format(read1,read2,sid,cpu_job,ram_job))

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print("Reads cannot be in multiple locations. Exiting.")
            sys.exit()
        else:
            pass

    #start multiprocessing
    pool = mp.Pool(processes=jobs)
    print("Begining assembly of reads:\n Number of Jobs: {0}\n CPUs/Job: {1}".format(jobs,cpu_job))

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Assembly\n')
        #begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/shovill:1.0.4',cmd,'/data',{read_path:"/data",os.path.join(outdir,'assemblies'):"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished Assembling Reads")



#annotation function
def annotate_assemblies(jobs,cpu_job,outdir):
    logfile = os.path.join(outdir,'annotation.log')
    input_path = os.path.join(outdir,"assemblies")
    #annotation
    annotated_path = os.path.join(outdir,"annotated")
    checkexists(annotated_path)

    #get assemblies
    fastas = list(getfiles(input_path))[0]
    #setup command list for annotating
    cmds = []
    for path in fastas:
        #main command
        assembly_file = os.path.basename(path)
        assembly_dir = os.path.basename(os.path.dirname(path))
        sid = os.path.basename(path).split('.')[0]
        cmds.append('prokka --compliant --outdir /output/{0} --prefix {1} --cpus {2} {3}/{4}'.format(sid,sid,cpu_job,assembly_dir,assembly_file))

    #start multiprocessing
    pool = mp.Pool(processes=jobs)
    print("Begining annotation of assemblies:\n Number of Jobs: {0}\n CPUs/Job: {1}".format(jobs,cpu_job))

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Annotating\n')
        #begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/prokka:1.13',cmd,'/data',{input_path:"/data",os.path.join(outdir,'annotated'):"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')

    #move finished .gff up a dir
    for root,dirs,files in os.walk(annotated_path):
        for file in files:
            if '.gff' in file:
                current_path = os.path.join(root,file)
                destination_path = os.path.join(os.path.dirname(root),file)
                os.rename(current_path,destination_path)

    print("Finished Annotating Assemblies")

#align assemblies function
def align(jobs,cpu_job,outdir,method='mafft'):
    cpus = jobs * cpu_job
    logfile = os.path.join(outdir,'alignment.log')
    input_path = os.path.join(outdir,"annotated")
    #remove alignment dir if it exists
    shutil.rmtree(os.path.join(outdir,'alignment'),ignore_errors=True)
    #alignment
    command = "sh -c 'roary -e -{0} -p {1} -f /output/alignment *.gff'".format(method,cpus)
    print("Begining alignment of core genes:\n Number CPUs: {0}".format(cpus))

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Aligning Core Genes\n')
        stdout = cd.call('staphb/roary:3.12.0',command,'/data',{input_path:"/data",outdir:"/output"})
        outlog.write('-----------\n')
        outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished Annotating Assemblies")

#create tree
def build_tree(outdir,model='GTR+G'):
    logfile = os.path.join(outdir,'tree.log')
    input_path = os.path.join(outdir,"alignment")

    #remove tree dir if it exists
    shutil.rmtree(os.path.join(outdir,'cg_tree'),ignore_errors=True)
    #create path
    cg_path = os.path.join(outdir,"cg_tree")
    checkexists(cg_path)
    shutil.copyfile(os.path.join(input_path,'core_gene_alignment.aln'),os.path.join(cg_path,'core_gene_alignment.aln'))

    #iqtree command
    command = "sh -c 'iqtree -s core_gene_alignment.aln -m {0} -bb 1000 '".format(model)
    print("Beginning to buld the Phylogeny")

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Building Tree\n')
        stdout = cd.call('staphb/iqtree:1.6.7',command,'/data',{cg_path:"/data"})
        outlog.write('-----------\n')
        outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished Bulding Phylogeny")

# ------------------------------------------------------

def core_genome(jobs,cpu_job,outdir):

    status,na = check_update_status(outdir,'','cg')
    if status == 'start':
        print("Starting the core-genome process.")
        print("Assembling reads using Shovill.")
        assemble_reads(jobs,cpu_job,outdir)
        check_update_status(outdir,"assemble",'cg')
        status = "assemble"

    if status == 'assemble':
        print("Annotating assemblies using Prokka.")
        annotate_assemblies(jobs,cpu_job,outdir)
        check_update_status(outdir,'annotate','cg')
        status = 'annotate'

    if status == 'annotate':
        print("Aligning Core Gene Set")
        align(jobs,cpu_job,outdir)
        check_update_status(outdir,'5','cg')
        status = 'align'

    if status == 'align':
        print("Building the Phylogeny")
        build_tree(outdir)
        check_update_status(outdir,'tree','cg')
        in_path = [outdir,'cg_tree','core_gene_alignment.aln.contree']
        out_path = [outdir,'core_genome_tree.tree']
        shutil.copyfile(os.path.join(*in_path),os.path.join(*out_path))
        status = 'tree'
    check_update_status(outdir,'done')
