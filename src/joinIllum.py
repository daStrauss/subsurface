'''
Created on Oct 3, 2012

@author: dstrauss
'''


import numpy as np
import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]

D = {'solverType':'splitField', 'flavor':'both', 'numRuns':1, 'expt':'joinEm', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    
    noFreqs = np.array(8) 
    bkg = 0

    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs))
    D['inc'] = np.array([45*np.pi/180.0])
    D['bkgNo'] = bkg+100;
    D['maxIter'] = 10
            
    return D