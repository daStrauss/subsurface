'''
Created on Oct 25, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'projection', 'flavor':'TE', 'numRuns':1, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkno = np.meshgrid(np.logspace(-8,0,20), np.arange(5))
    rhos = rhos.flatten()
    bkno = bkno.flatten()
    
    
    
    D['rho'] = rhos[parseNumber%100] 
    # D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = bkno[parseNumber] + 100
    
    D['freqs'] = np.array([1e3, 3e3, 13e3, 25e3])  
    D['inc'] = np.array([75])*np.pi/180
    D['numProcs'] = 4
    D['maxIter'] = 2
    
    return D