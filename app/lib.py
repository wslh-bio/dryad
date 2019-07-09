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

def check_update_status(path,status=''):
    path = os.path.abspath(path)
    if status != '':
        status_file = os.path.join(path,'status')
        with open(status_file,'a') as outstat:
            outstat.write(status+'\n')
    else:
        if "dryad-" in os.path.basename(path):
            status_file = os.path.join(path,"status")
            with open(status_file,'r') as instat:
                code = 0
                for line in instat:
                    code = int(line)
                if code != "done":
                    return code, path
        for root,dirs,files in os.walk(path):
            for dir in dirs:
                if "dryad-" in dir:
                    dryad_path = os.path.join(root,dir)
                    status_file = os.path.join(dryad_path,"status")
                    with open(status_file,'r') as instat:
                        code = 0
                        for line in instat:
                            code = int(line)
                        if code != "done":
                            return code, dryad_path
        return 0,""



    if not os.path.isdir(path):
        os.mkdir(path)
        return False
    else:
        return True

def getfiles(path):
    for root,dirs,files in os.walk(path):
        #scan path and look for files

        fastq_files = []
        bam_files = []
        gff_files = []

        for file in files:
            if "fastq.gz" in file:
                fastq_files.append(os.path.join(root,file))
            if "bam" in file:
                bam_files.append(os.path.join(root,file))
            if "gff" in file:
                gff_files.append(os.path.join(root,file))

        if len(fastq_files) > 0:
            fastq_files.sort()
            if len(fastq_files) % 2 != 0:
                print('There is an uneven number of read pairs in {0}. Exiting.'.format(path))
                sys.exit()
            paired_reads = []
            [paired_reads.append([x,y]) for x,y in zip(fastq_files[0::2],fastq_files[1::2])]
            yield paired_reads

        yield bam_files
        yield gff_files
