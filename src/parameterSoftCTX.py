'''
Created on Oct 18, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastSoftX', 'flavor':'TE', 'numRuns':500, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis = np.meshgrid(np.logspace(0,4,10), np.logspace(-8,-2,10))
    rhos = rhos.flatten()
    xis = xis.flatten()
    
    
    
    D['rho'] = rhos[parseNumber%100] 
    D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = int(parseNumber/100) + 100
    
    D['freqs'] = np.array([1e3, 3e3, 13e3, 25e3])  
    D['inc'] = np.array([75])*np.pi/180
    D['numProcs'] = 4
    
    return D