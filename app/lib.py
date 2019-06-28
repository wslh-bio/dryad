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


def getfiles(path):
    for root,dirs,files in os.walk(path):
        #scan path and look for files

        fastq_files = []
        bam_files = []

        for file in files:
            if "fastq.gz" in file:
                fastq_files.append(file)
            if "bam" in file:
                bam_files.append(file)

        if len(fastq_files) > 0:
            fastq_files.sort()
            if len(fastq_files) % 2 != 0:
                print('There is an uneven number of read pairs in {0}. Exiting.'.format(path))
                sys.exit()
            paired_reads = []
            [paired_reads.append([x,y]) for x,y in zip(fastq_files[0::2],fastq_files[1::2])]
            yield paired_reads

        yield bam_files
