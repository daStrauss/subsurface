'''
Created on Jun 18, 2012

@author: dstrauss
'''

import xml.etree.ElementTree as xml
import subprocess
import time
import sys
import os
import doFolders
    
def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return '''
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


def submitJob(cmd):
    pipeOut = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    f = pipeOut.stdout.read()
    f = f[:-1]
    print f + ' ' + time.asctime(time.localtime())
    time.sleep(1)
    return f



def main():
    if len(sys.argv) == 3:
        startIx = int(sys.argv[2])
    else:
        startIx = 0
    
    prSpec = __import__(sys.argv[1])
    finalIx = prSpec.D['numRuns']
    
    
    for ix in range(startIx,finalIx):
        doFolders.ensureFolders(prSpec.D, ix)
        jobTitle = 'run' + sys.argv[1] + repr(ix)
        fileName = 'sub' + sys.argv[1] + repr(ix) + '.pbs'
        
        fid = open(fileName, 'w')
        fid.write('mpiexec -wdir /shared/users/dstrauss/subsurface/src python coordinate.py ' + sys.argv[1] + ' ' + repr(ix))
        fid.close()
        cmd = ['qsub', '-N', jobTitle, '-l' , 'walltime=10:00:00', '-l','nodes=2:ppn=8', '-l', 'nice=0', fileName]
        
        print cmd
        ppid = submitJob(cmd)
        waitForExit(ppid)
        
if __name__=='__main__':
    main()
