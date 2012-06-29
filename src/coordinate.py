'''
Created on Jun 6, 2012

@author: dstrauss

Putting together a script to just simply loop through all of these materials and run on the cluster.
'''


import solveADMM
import sys
import os

# trialNs = 100

def main():
    ''' simple main routine '''
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        ix = int(sys.argv[2])
        outDir = 'splitField/' + sys.argv[3] + '/trial' + repr(ix) + '/'
        if sys.argv[3] == 'TE':
            solveADMM.semiParallel('splitField', sys.argv[3], rho=1500, xi =2e-3, \
                               uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1), outDir=outDir)
        elif sys.argv[3] == 'TM':
            solveADMM.semiParallel('splitField', sys.argv[3], rho=0.019307, xi =1.3895e-3, \
                               uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1), outDir=outDir)
            
    elif sys.argv[1] == 'contrastX':
        ix = int(sys.argv[2])
        outDir = 'contrastX/trial' + repr(ix)  + '/'
        
        solveADMM.semiParallel('contrastX', sys.argv[3], rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir = outDir)

    elif sys.argv[1] == 'sba':
        ix = int(sys.argv[2])
        outDir = 'sba/trial'+repr(ix) + '/'
        solveADMM.semiParallel('sba', sys.argv[3], rho=0.005, xi=0.9, uBound=0.05, lmb=0, bkgNo=(ix+1), outDir=outDir)
        
    elif sys.argv[1] == 'biconvex':
        ix = int(sys.argv[2])
        outDir = 'biconvex/trial' + repr(ix) + '/'
        solveADMM.semiParallel('biconvex', sys.argv[3], rho=0.001, xi=1e-5, lmb=0, uBound=0.05,bkgNo=(ix+1), outDir=outDir)
        
    else: 
        print 'I think you asked for the wrong thing:'
 
if __name__ == "__main__":
    main()
