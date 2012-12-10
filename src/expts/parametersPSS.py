'''
Created on Nov 28, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'phaseSplitSingle', 'flavor':'TE', 'numRuns':960, 'expt':'paramPSS', 'numProcs':1}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis, fr, bkgLocal = np.meshgrid(np.logspace(-4,4,8), np.linspace(0,1,8), np.array([1e3, 1e4, 1e5]), range(5))
    rhos = rhos.flatten()
    xis = xis.flatten()
    bkgLocal = bkgLocal.flatten()
    fr = fr.flatten()
        
    D['freqs'] = np.array(fr[parseNumber])
    D['numProcs'] = 1
    D['numSensors'] = 3100

    D['lam'] = 0
#    D['rho'] = 1e-3    
    D['inc'] = np.array([75*np.pi/180])
    D['maxIter'] = 200
    
    D['rho'] = rhos[parseNumber] 
    D['xi'] = xis[parseNumber]
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    return D