'''
Created on Jun 15, 2012

@author: dstrauss

Putting together a script to do a coordinated search and find the optimal parameters for the solver
(6/15) - Trying first with the semiParallel routine and the TM splitField formulation.

(6/21) -- adding the ability to separate each of the tasks into a separate job. Don't overload the queue.

'''


import solveADMM
import sys
import os
import numpy as np
import subprocess

def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return '''
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


def main():
    ''' simple main routine '''
    
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        assert os.path.exists('splitField')
        
        rhos, xis = np.meshgrid(np.logspace(-2,2,8), np.logspace(-3,0,8))
        rhos = rhos.flatten()
        xis = xis.flatten()
        
        for supIx in range(trialNs):
            for ix, r in enumerate(rhos):
                outDir = 'splitField/mat' + repr(supIx) + '/trial' + repr(ix) + '/'
                assert os.path.exists(outDir)            
                solveADMM.semiParallel('splitField', 'TM', rho=r, xi=xis[ix], \
                      uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1), outDir=outDir)
            
            
    elif sys.argv[1] == 'contrastX':
        if not os.path.exists('contrastX'):
            os.mkdir('contrastX')
            
        for ix in range(trialNs):
            outDir = 'contrastX/trial' + repr(ix)  + '/'
            if not os.path.exists(outDir):
                os.mkdir(outDir)

            solveADMM.semiParallel('contrastX', rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir = outDir)

    elif sys.argv[1] == 'sba':
        if not os.path.exists('sba'):
            os.mkdir('sba')
        
        
        for ix in range(trialNs):
            outDir = 'sba/trial'+repr(ix) + '/'
            if not os.path.exists(outDir):
                os.mkdir(outDir)

            solveADMM.semiParallel('sba', rho=0.005, xi=0.9, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir=outDir)
    else: 
        print 'I think you asked for the wrong thing:'
 
if __name__ == "__main__":
    main()
