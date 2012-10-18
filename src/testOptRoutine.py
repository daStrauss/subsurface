'''
Created on Oct 17, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'contrastSoftX', 'flavor':'TE', 'numRuns':1, 'expt':'testSolver'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    D['rho'] = 0.001
    D['xi'] = 0.0001
    D['freqs'] = np.array([1e3])
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  100
    D['numProcs'] = 1
    D['maxIter'] = 5
    
    return D