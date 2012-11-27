'''
Created on Nov 14, 2012

@author: dstrauss

Simple routine to do optimization over 3D! 

'''


import numpy as np
D = {'solverType':'ctxdev', 'flavor':'TE3D', 'numRuns':3, 'expt':'huge', 'numProcs':1}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
    D['freqs'] = np.array([1e3])
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['numProcs'] = 1
    D['reg'] = 1e-8
    if parseNumber == 0:
        D['maxIter'] = 1
    elif parseNumber == 1:
        D['maxIter'] = 2
    else:
        D['maxIter'] = 10
        
        
    D['lmb'] = 1e-13
    
    return D