#author: Kelsey Florek
#email: kelsey.florek@slh.wisc.edu
#stripped down docker calling function

import docker
import os, sys
import time

def call(container,command,cwd='',paths={},remove=True):
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
        container_obj = client.containers.run(container,command,user=user,volumes=volumes,working_dir=cwd,remove=remove,detach=True)
    except:
        #loop through output as it is streamed
        for line in container_obj.logs(stream=True):
            output += line
    else:
        for line in container_obj.logs(stream=True):
            output += line
    #once container is finished return output as a string
    return output.decode('utf-8')
