'''
Created on Jun 21, 2012

@author: dstrauss

more notes! update in November:
I now know more about what I want this to do, so I'm going to do it.
I know how to parse the xml files to look at what jobs are currently running

Differences w.r.t. the maui/torque platform:
(1) returns from qsub are not valuable
(2) query all of qstat to get xml output
(3) search the jobnames that are running to see if they match
(4) return false otherwise.
(5) make the joblist just the names of the scripts
(6) jobname is really not necessary

'''

import xml.etree.ElementTree as xml
import subprocess
import time
import sys
import doFolders
import math
import importlib

class jobParm(object):
    def __init__(self, name, folder):
        self.name = name
        self.folder = folder
   
def checkForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return.
    Returns True if the job is still running with status R, returns false if job has completed
    or simply does not have status R. 
    update! jobname is now a class/structure as defined above!
    '''
    doExit = True
    pipeOut = subprocess.Popen(['qstat', '-xml'], stdout=subprocess.PIPE)
    f = pipeOut.stdout.read()
    # print f
    try:
        P = xml.fromstring(f)
        # print 'parsed string'
        # rtf = P.getroot() ? interesting it returns the root
        # print 'got root'
        ''' just have to check to see if it is in the joblist -- they disappear quickly'''
        for job in P.iter('job_list'):
            for z in job.iter('JB_name'):
                if z.text == jobName.name:
                    doExit = False

    except:
        print 'wrong string or something?'
    return doExit

  
    
def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return '''
    doExit = False
    while doExit == False:
        pipeOut = subprocess.Popen(['qstat', '-xml', jobName], stdout=subprocess.PIPE)
        f = pipeOut.stdout.read()
        P = xml.fromstring(f)
        for z in P.getiterator():
            if z.tag == 'job_state':
                if z.text == 'R':
                    pass
                    # print 'Still Running ' + jobName
                elif z.text == 'C':
                    doExit = True
                    print 'Finished ' + jobName
        
        time.sleep(5)


def submitJob(cmd):
    pipeOut = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    f = pipeOut.stdout.read()
    f = f[:-1]
    print f + ' ' + time.asctime(time.localtime())
    # time.sleep(1)
    return f

def main():
    ''' main routine: creates a runList and a jobList. JobList are currently running, 
    runList are too run
    '''
    ''' parse the starting index, should be the second argument '''
    if len(sys.argv) >= 3:
        startIx = int(sys.argv[2])
    else:
        startIx = 0
    
    ''' get the number of simultaneous workers, i.e. pool size '''
    if len(sys.argv) == 4:
        numWorkers = int(sys.argv[3])
    else:
        numWorkers = 1
        
    ''' import the management script '''
    prSpec = importlib.import_module('expts.' + sys.argv[1])
    finalIx = prSpec.D['numRuns']
    
    ''' creat a list of jobs to run '''
    runList = range(startIx,finalIx)
    
    ''' create a pool of jobs '''
    jobList = list()
    
    while len(runList) > 0:
        
        ''' populate the joblist '''
        if len(jobList) < numWorkers:
            # launch worker
            ix = runList.pop(0)
            mvFolder = doFolders.ensureFolders(prSpec.D, ix)
            
            lclD = prSpec.getMyVars(ix, prSpec.D)
            nProcs = lclD['numProcs']
            nNodes = int(math.ceil(nProcs/8.0))
            
            ''' create a submit.pbs file '''
            fileName = 'sub' + sys.argv[1] + repr(ix) + '.pbs'
        
            fid = open(fileName, 'w')
            fid.write('mpiexec -wdir /home/dstrauss/subsurface/src python coordinate.py ' + sys.argv[1] + ' ' + repr(ix))
            fid.close()
            
            cmd = ['qsub', '-pe', 'orte', repr(nProcs), '-cwd', fileName]    
            print cmd
                
            submitJob(cmd)
            
            nJob = jobParm(fileName,mvFolder)
            jobList.append(nJob)
            
        else:
            time.sleep(5)
    
        for jbs in jobList:
            if checkForExit(jbs):
                cmd = ['mkdir', '-p', 'mnt/dstrauss/' + jbs.folder]
                print cmd
                submitJob(cmd)
                cmd = ['rsync', '-a', jbs.folder, '/mnt/dstrauss/' + jbs.folder]
                print cmd
                submitJob(cmd)
                cmd = ['rm', '-r', jbs.folder]
                print cmd
                submitJob(cmd)
                
                jobList.remove(jbs)
  
        
if __name__=='__main__':
    main()
