import os,sys
import multiprocessing as mp
import psutil

from lib import getfiles,checkexists

def core_genome(input_path,jobs,cpu_job,outdir):

    print("Starting the core-genome process.")
    logfile = os.path.join(outdir,'core_genome.log')

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
        sid = os.path.basename(read_pair[0]).split('_')[0]
        cmds.append('shovill --R1 {0} --R2 {1} --outdir {2} --cpus {3} --ram {4}'.format(read_pair[0],read_pair[1],sid,cpu_job,ram_job))

        if read_path == '':
            read_path = os.path.dirname(os.path.abspath(read_pair[1]))
        elif read_path != os.path.dirname(os.path.abspath(read_pair[1])):
            print("Reads cannot be in multiple locations. Exiting.")
            sys.exit()
        else:
            pass

    #start multiprocessing
    pool = mp.Pool(processes=jobs)
    print("Begining assembly of reads:\n Number of Jobs: {0}\n CPUs/Job: {1}".format(jobs,cpu))

    #denote logs
    with open(logfile,'a') as outlog:
        outlog.write('***********\n')
        outlog.write('Assembly\n')
        #begin multiprocessing
        results = pool.starmap_async(cd.call,[['staphb/trimmomatic:0.39',cmd,'/data',{read_path:"/data",os.path.join(outdir,'trimmed'):"/output"}] for cmd in cmds])
        stdouts = results.get()
        for stdout in stdouts:
            outlog.write('-----------\n')
            outlog.write(stdout)
        #denote end of logs
        outlog.write('***********\n')
    print("Finished Assembling Reads")

    #annotate
    #align
    #constrct_tree
