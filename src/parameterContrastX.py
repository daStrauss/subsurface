'''
Created on Jul 30, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TM', 'numRuns':1, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(5))
    rhos = rhos.flatten()
    bkgLocal = bkgLocal.flatten()
    
    
    
    D['rho'] = rhos[parseNumber] 
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 16
    return D