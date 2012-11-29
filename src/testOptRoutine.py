'''
Created on Oct 17, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'phaseSplitSingle', 'flavor':'TE', 'numRuns':1, 'expt':'testSolver'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
    D['freqs'] = np.array([1e3])
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  100
    D['numProcs'] = 1
    D['maxIter'] = 100
    D['numSensors'] = 3010
    D['rho'] = 100.0
    
    return D