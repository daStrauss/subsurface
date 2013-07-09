'''
Created on Nov 8, 2012
Copyright © 2013
The Board of Trustees of The Leland Stanford Junior University.
All Rights Reserved

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
@author: dstrauss

This ought to be identical to manageTasks, except that it
includes a flag in the MPI function calls to work properly on the new
nodes of our cluster. 
'''

import xml.etree.ElementTree as xml
import subprocess
import time
import sys
import doFolders
import math
import importlib
    
def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return.
    This program takes full control and does not return until the job has completed. '''
    doExit = False
    while doExit == False:
        pipeOut = subprocess.Popen(['qstat', '-x', jobName], stdout=subprocess.PIPE)
        f = pipeOut.stdout.read()
        try:
            P = xml.fromstring(f)
            for z in P.getiterator():
                if z.tag == 'job_state':
                    if z.text == 'R':
                        pass
                        # print 'Still Running ' + jobName
                    elif z.text == 'C':
                        doExit = True
                        print 'Finished ' + jobName
        except:
            print 'wrong string or something?'
        
        time.sleep(5)
        
        
def checkForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return.
    Returns True if the job is still running with status R, returns false if job has completed
    or simply does not have status R. '''
    doExit = False
    pipeOut = subprocess.Popen(['qstat', '-x', jobName], stdout=subprocess.PIPE)
    f = pipeOut.stdout.read()
    try:
        P = xml.fromstring(f)
        for z in P.getiterator():
            if z.tag == 'job_state':
                if z.text == 'R':
                    pass
                        # print 'Still Running ' + jobName
                elif z.text == 'C':
                    doExit = True
                    print 'Finished ' + jobName
    except:
        print 'wrong string or something?'
    return doExit


def submitJob(cmd):
    ''' Basic routine to submit/execute a line of code given in cmd and return whatever return
    values are spit out by the system.'''
    pipeOut = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    f = pipeOut.stdout.read()
    f = f[:-1]
    print f + ' ' + time.asctime(time.localtime())
    time.sleep(1)
    return f



def main():
    ''' main routine to run one of my python submission routines. It takes up to three arguments
    in the command line. It wants an iteration control script, a start index, and a number of submissions
    to maintain. From there, it just goes'''
    # problem specific issues
    if len(sys.argv) >= 3:
        startIx = int(sys.argv[2])
    else:
        startIx = 0
    
    if len(sys.argv) == 4:
        numWorkers = int(sys.argv[3])
    else:
        numWorkers = 1
    

    # import my runtime iteration control script
    prSpec = importlib.import_module('expts.' + sys.argv[1])
    finalIx = prSpec.D['numRuns']
    
    # queue of indexes still to run
    runList = range(startIx,finalIx)
    # list of open/running jobs
    # jobList = list(range(numWorkers))
    jobList = list()
    
    # main loop: wait for all potential jobs to be executed
    while len(runList) > 0:
        # populate the job list
        if len(jobList) < numWorkers:
            # launch worker
            ix = runList.pop(0)
            doFolders.ensureFolders(prSpec.D, ix)
            
            lclD = prSpec.getMyVars(ix, prSpec.D)
            nProcs = lclD['numProcs']
            # nNodes = int(math.ceil(nProcs/8.0))
            # nProcs = 6
            
            if nProcs > 0:
                jobTitle = 'run' + sys.argv[1] + repr(ix)
                fileName = 'sub' + sys.argv[1] + repr(ix) + '.pbs'
        
                fid = open(fileName, 'w')
                fid.write('mpiexec -n ' + repr(nProcs) + ' --mca btl_tcp_if_include "p3p2" -wdir /shared/users/dstrauss/subsurface/src python coordinate.py ' + sys.argv[1] + ' ' + repr(ix))
                fid.close()
                cmd = ['qsub', '-N', jobTitle, '-l' , 'walltime=20:00:00', '-l','nodes=' + repr(nProcs), '-l', 'nice=0', '-q','batchnew', fileName]        
                print cmd
            
                jobList.append(submitJob(cmd))
            else:
                print 'skipping ' + repr(ix)
                
        else:
            time.sleep(5)   
         
         # check for completion   
        for jbs in jobList:
            if checkForExit(jbs):
                jobList.remove(jbs)
                
        
if __name__=='__main__':
    main()
