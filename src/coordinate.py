'''
Created on Jun 6, 2012

@author: dstrauss

Putting together a script to just simply loop through all of these materials and run on the cluster.
'''


import solveADMM
import sys
import os

trialNs = 5

def main():
    ''' simple main routine '''
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        if not os.path.exists('splitField'):
            os.mkdir('splitField')

        for ix in range(10,15):
            outDir = 'splitField/trial' + repr(ix) + '/'
            if not os.path.exists(outDir):
                os.mkdir(outDir)
            
            solveADMM.semiParallel('splitField', rho=1500, xi =2e-4, \
                                   uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1), outDir =outDir)
            
            
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
