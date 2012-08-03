'''
Created on Aug 3, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':2000, 'expt':'bkgSig'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(100))
    rhos = rhos.flatten()
    bkgLocal = bkgLocal.flatten()
    
    
    
    D['bkgSig'] = rhos[parseNumber] 
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 16
    return D