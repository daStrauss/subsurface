'''
Created on Jun 6, 2012

@author: dstrauss

Putting together a script to just simply loop through all of these materials and run on the cluster.
'''


import solveADMM
import sys
import os

trialNs = 100

def main():
    ''' simple main routine '''
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        ix = int(sys.argv[2])
        outDir = 'splitField/trial' + repr(ix) + '/'
        solveADMM.semiParallel('splitField', 'TE', rho=1500, xi =2e-3, \
                               uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1), outDir=outDir)
            
    elif sys.argv[1] == 'contrastX':
        ix = int(sys.argv[2])
        outDir = 'contrastX/trial' + repr(ix)  + '/'
        
        solveADMM.semiParallel('contrastX', 'TE', rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir = outDir)

    elif sys.argv[1] == 'sba':
        ix = int(sys.argv[2])
        outDir = 'sba/trial'+repr(ix) + '/'
        solveADMM.semiParallel('sba', 'TE', rho=0.005, xi=0.9, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir=outDir)
        
    else: 
        print 'I think you asked for the wrong thing:'
 
if __name__ == "__main__":
    main()
