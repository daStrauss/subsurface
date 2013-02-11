'''
Created on Feb 11, 2013

@author: dstrauss
'''

import numpy as np


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':800, 'expt':'paramPHS', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis, bkgLocal = np.meshgrid(np.logspace(-4,4,10), np.linspace(-13,-5,10), range(8))
    rhos = rhos.flatten()
    xis = xis.flatten()
    bkgLocal = bkgLocal.flatten()
        
    D['lam'] = 0
#    D['rho'] = 1e-3    

    D['rho'] = rhos[parseNumber] 
    D['xi'] = xis[parseNumber]
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    return D