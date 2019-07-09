import os,sys
import multiprocessing as mp
import psutil

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
    fastqs,bam = getfiles(input_path)
    cmds = []
    read_path = ''
    for read_pair in fastqs:
        #main command
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])
        sid = os.path.basename(read_pair[0]).split('_')[0]
        #TODO get renaming of reads working
        cmds.append('shovill --R1 /data/{0} --R2 /data/{1} --outdir /output/{2} --cpus {3} --ram {4} && mv /output/{2}/contigs.fa /output/{2}/{2}.fa '.format(read1,read2,sid,cpu_job,ram_job))

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
    input_path = os.path.join(outdir,'assemblies')

    #annotation
    annotated_path = os.path.join(outdir,"annotated")
    checkexists(annotated_path)

    #get trimmed reads
    #TODO test GFF searching
    fastqs,bam = getfiles(input_path)
    cmds = []
    read_path = ''
    for read_pair in fastqs:
        #main command
        read1 = os.path.basename(read_pair[0])
        read2 = os.path.basename(read_pair[1])
        sid = os.path.basename(read_pair[0]).split('_')[0]
        cmds.append('shovill --R1 /data/{0} --R2 /data/{1} --outdir /output/{2} --cpus {3} --ram {4}'.format(read1,read2,sid,cpu_job,ram_job))

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

#align assemblies function
def align():
    pass

#create tree
def build_tree():
    pass

# ------------------------------------------------------

def core_genome(jobs,cpu_job,outdir):

    status,na = check_update_status(outdir)

    if status == 2:
        print("Starting the core-genome process.")
        print("Assembling reads using Shovill.")
        assemble_reads(jobs,cpu_job,outdir)
        check_update_status(outdir,"3")
        status = 3

    if status == 3:
        print("Annotating assemblies using Prokka.")
