#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#miscellaneous functions

import os

#returns (num jobs,num cpus per job)
def cpu_count(num):
    if num <=1:
        return 1,1
    elif num <= 5:
        return 1,num
    else:
        results = []
        for n in range(2,num):
            if num % n == 0:
                results.append(n)
        if len(results) < 2:
            return cpu_count(num-1)
        index = int(len(results)/2)
        return results[index],int(num/results[index])

def checkexists(path):
    path = os.path.abspath(path)
    if not os.path.isdir(path):
        os.mkdir(path)
        return False
    else:
        return True

def check_update_status(path,status='',pipe=''):
    if pipe:
        path += '_'+pipe
    path = os.path.abspath(path)
    if status != '':
        status_file = os.path.join(path,'status')
        with open(status_file,'a') as outstat:
            outstat.write(status+'\n')
    else:
        if "dryad-" in os.path.basename(path):
            status_file = os.path.join(path,"status")
            with open(status_file,'r') as instat:
                code = ''
                for line in instat:
                    code = line
                if code != "done":
                    return code, path
        for root,dirs,files in os.walk(path):
            for dir in dirs:
                if "dryad-" in dir:
                    dryad_path = os.path.join(root,dir)
                    status_file = os.path.join(dryad_path,"status")
                    with open(status_file,'r') as instat:
                        code = ''
                        for line in instat:
                            code = line
                        if code != "done":
                            return code, dryad_path
        return '',""



    if not os.path.isdir(path):
        os.mkdir(path)
        return False
    else:
        return True

def getfiles(path):
    #scan path and look for files
    fastq_files = []
    fasta_files = []
    gff_files = []
    for root,dirs,files in os.walk(path):
        for file in files:
            if ".fastq.gz" in file:
                fastq_files.append(os.path.join(root,file))
            if ".fa" in file and "spades" not in file:
                fasta_files.append(os.path.join(root,file))
            if ".gff" in file:
                gff_files.append(os.path.join(root,file))

    if len(fastq_files) > 0:
        fastq_files.sort()
        if len(fastq_files) % 2 != 0:
            print('There is an uneven number of read pairs in {0}. Exiting.'.format(path))
            sys.exit()
        paired_reads = []
        [paired_reads.append([x,y]) for x,y in zip(fastq_files[0::2],fastq_files[1::2])]
        yield paired_reads

    if len(fasta_files) > 0:
        yield fasta_files

    if len(gff_files) > 0:
        yield gff_files
