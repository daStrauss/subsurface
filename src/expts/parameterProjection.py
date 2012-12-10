'''
Created on Oct 25, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'projection', 'flavor':'TE', 'numRuns':10, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos = np.logspace(-8,0,10)
    
#     , bkno = np.meshgrid(np.logspace(-8,0,20), np.arange(5))
    # rhos = rhos.flatten()
    # bkno = bkno.flatten()
    
    
    
    D['rho'] = rhos[parseNumber] 
    # D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = 100
    
    D['freqs'] = np.array([1e3, 3e3, 13e3, 25e3])  
    D['inc'] = np.array([75])*np.pi/180
    D['numProcs'] = 4
    D['maxIter'] = 100

    
    return D