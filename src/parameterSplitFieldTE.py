'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':500, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis = np.meshgrid(np.logspace(2,4,10), np.logspace(-4,-2,10))
    rhos = rhos.flatten()
    xis = xis.flatten()
    
    
    
    D['rho'] = rhos[parseNumber%100] 
    D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = int(parseNumber/100) + 100
    
    return D