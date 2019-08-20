#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#stripped down docker calling function

import docker
import os, sys
import time
import signal

def call(container,command,cwd='',paths={},remove=True,cpu_set='',sig_default=True):
    ###set signal handler to default if we are part of a subprocess
    if sig_default:
        signal.signal(signal.SIGINT,signal.SIG_DFL)
    ###access docker environment
    client = docker.from_env()

    ###get the effective user and group id's
    user = str(os.geteuid())+':'+str(os.getegid())

    ###setup mount point paths
    #{"/path/outside":"/path/incontainer"}
    volumes = {}
    if paths:
        for key in paths.keys():
            volumes[key] = {'bind':paths[key],'mode':'rw'}

    ###run the container
    #create empty variable for holding byte object output for the container logs
    output = b''
    #try block to run the container
    try:
        ###format cpu set for control of cpus
        if cpu_set:
            cpu_set = int(str(int(cpu_set))+"00000")
            container_obj = client.containers.run(container,command,user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True,cpu_period=100000,cpu_quota=cpu_set,labels={"prog":"dryad"})
        else:
            container_obj = client.containers.run(container,command,user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True,labels={"prog":"dryad"})
    except:
        #loop through output as it is streamed
        for line in container_obj.logs(stream=True):
            output += line
    else:
        for line in container_obj.logs(stream=True):
            output += line
    #once container is finished return output as a string
    return output.decode('utf-8')
