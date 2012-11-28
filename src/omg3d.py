'''
Created on Nov 14, 2012

@author: dstrauss

Simple routine to do optimization over 3D! 

'''


import numpy as np
D = {'solverType':'ctxdev', 'flavor':'TE3D', 'numRuns':4, 'expt':'huge', 'numProcs':1}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
    D['freqs'] = np.array([1e3])
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['numProcs'] = 1
    D['reg'] = 1e-6
    if parseNumber == 0:
        D['maxIter'] = 1
    elif parseNumber == 1:
        D['maxIter'] = 2
    elif parseNumber == 2:
        D['maxIter'] = 10
    elif parseNumber == 3:
        D['maxIter'] = 10
        D['reg'] = 0.0   
        
    elif parseNumber == 4:
        D['bkgSig'] = 1e-5
        D['maxIter'] = 30
        D['freqs'] = np.array([1e5])
        D['reg'] = 0.0 
        D['rho'] = 100.0
        
        
    D['lmb'] = 1e-5
    
    return D