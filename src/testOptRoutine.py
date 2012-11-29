'''
Created on Oct 17, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':2, 'expt':'testSolver'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
#    D['freqs'] = np.array([1e3])
#    D['inc'] = np.array([75*np.pi/180])
#    D['bkgNo'] =  100
#    D['numProcs'] = 1
#    D['maxIter'] = 100
#    D['numSensors'] = 3010
#    D['rho'] = 0.10
    if parseNumber == 0:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([45*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 1e-3
        D['xi'] = 1e-12
    else:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([45*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 5e-4
        D['xi'] = 1e-12
    
    return D