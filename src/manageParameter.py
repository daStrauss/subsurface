'''
Created on Jun 21, 2012

@author: dstrauss
'''

import subprocess
import xml.etree.ElementTree as xml
import time
import os
import sys


def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return 
        This only works well on nansen'''
    doExit = False
    while doExit == False:
        pipeOut = subprocess.Popen(['qstat', '-x', jobName], stdout=subprocess.PIPE)
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
    time.sleep(1)
    return f

def ensureFolders(baseTag,flavor):
    ''' here's where the folder making mechanism goes! '''
    
    if not os.path.exists(baseTag):
        os.mkdir(baseTag)
    
    #mid = baseTag + '/mat' + repr(supIx) + '/'
    mid = baseTag + '/' + flavor + '/'
    if not os.path.exists(mid):
        os.mkdir(mid)
        
    for ix in range(100):
        outDir = mid + '/prmTrial' + repr(ix) + '/'
        print outDir
        if not os.path.exists(outDir):
            os.mkdir(outDir)
            
        figDir = outDir + 'Figs/'
        datDir = outDir + 'Data/'
        if not os.path.exists(figDir):
            os.mkdir(figDir)
        if not os.path.exists(datDir):
            os.mkdir(datDir)


def main():
    ''' wrap around for making scripts and submiting them and waiting '''
    for ix in range(0,64):
        jobTitle = 'paramS' + sys.argv[1] + repr(ix)
        fileName = 'paramS' + sys.argv[1] + '.pbs'
        
        if sys.argv[1] == 'sba':
            fid = open(fileName, 'w')
            fid.write('mpiexec -npernode 8 -wdir /shared/users/dstrauss/subsurface/src python coordinate.py ' + sys.argv[1] + ' ' + repr(ix) + ' TE')
            fid.close()
            cmd = ['qsub', '-N', jobTitle, '-l' , 'walltime=10:00:00', '-l','nodes=8:ppn=8', fileName]
        
        elif sys.argv[1] == 'splitField':           
            fid = open(fileName, 'w')
            fid.write('mpiexec -npernode 8 -wdir /shared/users/dstrauss/subsurface/src python parameterSearch.py ' + sys.argv[1] + ' ' + repr(ix)) 
            fid.close()
            cmd = ['qsub', '-N', jobTitle, '-l' , 'walltime=10:00:00', '-l','nodes=2:ppn=8', fileName]
        print cmd
        ppid = submitJob(cmd)
        waitForExit(ppid)
