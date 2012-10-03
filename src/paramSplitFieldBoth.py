'''
Created on Oct 3, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'both', 'numRuns':500, 'expt':'intParameters', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis = np.meshgrid(np.logspace(2,4,10), np.logspace(-4,-2,10))
    rhos = rhos.flatten()
    xis = xis.flatten()
    
    noFreqs = np.array(8) 
    bkg = 0

    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs))
    D['inc'] = np.array([45*np.pi/180.0])
    
    D['rho'] = rhos[parseNumber%100] 
    D['xi'] = xis[parseNumber%100]
    D['bkgNo'] = int(parseNumber/100) + 100
    
    return D