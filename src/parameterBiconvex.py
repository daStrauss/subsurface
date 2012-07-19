'''
Created on Jul 17, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'biconvex', 'flavor':'TE', 'numRuns':100, 'expt':'intParameters'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkgLocal = np.meshgrid(np.logspace(1,5,20), range(5))
    rhos = rhos.flatten()
    bkgLocal = bkgLocal.flatten()
    
    
    
    D['rho'] = rhos[parseNumber] 
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    
    return D