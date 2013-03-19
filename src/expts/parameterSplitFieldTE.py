'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TM', 'numRuns':800, 'expt':'intParameters', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis = np.meshgrid(np.logspace(1,5,10), np.logspace(-5,-1,10))
    rhos = rhos.flatten()
    xis = xis.flatten()
    
    
    
    D['rho'] = rhos[parseNumber%100] 
    D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = int(parseNumber/100) + 100
    
    return D