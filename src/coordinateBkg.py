'''
Created on Jun 30, 2012

@author: dstrauss
Want to coordinate a search through 10 background profiles and 10 sigmas.
Created on Jun 15, 2012

@author: dstrauss
'''
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

numXi = 10
numRho = 10
totStates = numXi*numRho


def getMyVars(parseNumber):
    '''routine to return the parameters to test at the current iteration.'''
    sigmaBkg, profileNo = np.meshgrid(np.logspace(-3,0,numRho), np.arange(10))
    sigmaBkg = sigmaBkg.flatten()
    profileNo = profileNo.flatten()
    return profileNo[parseNumber], sigmaBkg[parseNumber]
    
def main():
    ''' simple main routine '''
    
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    elif sys.argv[1] == 'splitField':
        assert os.path.exists('splitField')
        
        localProfileNo,localSigmaBkg = getMyVars(int(sys.argv[2]))
        
        outDir = 'splitField/' + 'TE/' + '/bkgTrial' + sys.argv[2] + '/'
        
        assert os.path.exists(outDir)            
        # solveADMM.semiParallel('splitField', 'TE', rho=140, xi=1e-3, \
        #       uBound = 0.05, lmb = 1e-8, bkgNo=localProfileNo, outDir=outDir, bkg=localSigmaBkg)
            

if __name__ == "__main__":
    main()
