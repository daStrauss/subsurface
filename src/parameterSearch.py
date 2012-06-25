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

numXi = 8
numRho = 8
totStates = numXi*numRho


def getMyVars(parseNumber):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis = np.meshgrid(np.logspace(-4,0,numRho), np.logspace(-5,-2,numXi))
    rhos = rhos.flatten()
    xis = xis.flatten()
    return rhos[parseNumber], xis[parseNumber]
    
def main():
    ''' simple main routine '''
    
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        assert os.path.exists('splitField')
        
        lRho,lXi = getMyVars(int(sys.argv[2]))
        
        outDir = 'splitField/' + 'TM/' + '/prmTrial' + sys.argv[2] + '/'
        
        assert os.path.exists(outDir)            
        solveADMM.semiParallel('splitField', 'TM', rho=lRho, xi=lXi, \
              uBound = 0.05, lmb = 1e-8, bkgNo=1, outDir=outDir)
            

if __name__ == "__main__":
    main()
