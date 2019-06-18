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
