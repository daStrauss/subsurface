'''
Created on Dec 7, 2012

@author: dstrauss


'''

import numpy as np


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':200, 'expt':'standard', 'numProcs':4}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber + 100
    D['rho'] = 0.5e-3
    D['xi'] = 1e-12
    D['numProcs'] = 4
    D['freqs'] = np.array([1e3, 3e3, 13e3, 50e3])  
    D['inc'] = np.array([75])*np.pi/180
    
    
     
    return D