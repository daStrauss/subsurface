'''
Created on Nov 8, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':245, 'expt':'newParams2', 'numProcs':2}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis, bkgLocal = np.meshgrid(np.logspace(-4,0,7), np.logspace(-16,-8,7), range(5))
    rhos = rhos.flatten()
    xis = xis.flatten()
    bkgLocal = bkgLocal.flatten()
    
        
    D['freqs'] = np.array([1e3, 25e3])
    D['numProcs'] = 2
    D['numSensors'] = 400

    D['lam'] = 0
#    D['rho'] = 1e-3    
    D['inc'] = np.array([75*np.pi/180])
    D['maxIter'] = 200

    
    
    D['rho'] = rhos[parseNumber] 
    D['xi'] = xis[parseNumber]
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 2
    return D