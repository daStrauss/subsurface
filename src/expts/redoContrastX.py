'''
Created on Jul 18, 2012

@author: dstrauss
'''

'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':613, 'expt':'standard', 'numProcs': 16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber
#    D['rho'] = 0.5e-3
#    D['xi'] = 1e-12
#    D['numProcs'] = 16
#    D['freqs'] = np.array([1e3, 3e3, 13e3, 50e3])  
#    D['inc'] = np.array([75])*np.pi/180
    
    
     
    return D