'''
Created on Nov 8, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':100, 'expt':'newParams6', 'numProcs':6}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(5))
    rhos = rhos.flatten()
    bkgLocal = bkgLocal.flatten()
    
    D['freqs'] = np.linspace(1e3,25e3,6)
    D['numProcs'] = 6
    D['numSensors'] = 400
        
    D['lam'] = 1e-8
    # D['rho'] = 1e-3    
    D['inc'] = np.array([75*np.pi/180])
    
    D['rho'] = rhos[parseNumber] 
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 6
    return D