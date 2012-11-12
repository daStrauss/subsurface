'''
Created on Nov 7, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'middleMan', 'flavor':'TE', 'numRuns':4, 'expt':'testThree'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    if (parseNumber == 0):
        D['freqs'] = np.array([1e3])
        D['numProcs'] = 1
        D['numSensors'] = 400
    elif (parseNumber == 1):
        D['freqs'] = np.array([1e3, 25e3])
        D['numProcs'] = 2
        D['numSensors'] = 400
    elif (parseNumber == 2):
        D['freqs'] = np.array([25e3])
        D['numProcs'] = 1
        D['numSensors'] = 400
    elif (parseNumber == 3):
        D['freqs'] = np.linspace(1e3,25e3,6)
        D['numProcs'] = 6
        D['numSensors'] = 400
        
    D['lam'] = 0.0
    D['rho'] = 0.001 
    D['xi'] = 1e-7 
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['maxIter'] = 20
    
    return D
