'''
Created on Feb 5, 2013

@author: dstrauss
'''


import numpy as np


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':2000, 'expt':'bkgSig'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(100))
    rhos = rhos.flatten()
    bkgLocal = bkgLocal.flatten()
    
    
    
    D['bkgSig'] = rhos[parseNumber] 
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['inc'] = np.array([75.0*np.pi/180.0])
    D['numProcs'] = 4
    D['rho'] = 1e-3
    D['xi'] = 1e-12
    
    return D