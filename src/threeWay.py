'''
Created on Nov 7, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':3, 'expt':'testThree'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    if (parseNumber == 0):
        D['freqs'] = np.array([1e3])
        D['numProcs'] = 1
    elif (parseNumber == 1):
        D['freqs'] = np.array([1e3, 25e3])
        D['numProcs'] = 2
    elif (parseNumber == 2):
        D['freqs'] = np.array([25e3])
        D['numProcs'] = 1
        
        
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['maxIter'] = 10
    
    return D
