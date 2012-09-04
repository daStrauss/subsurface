'''
Created on Sep 4, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':3, 'expt':'spfTest'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
#    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(100))
#    rhos = rhos.flatten()
#    bkgLocal = bkgLocal.flatten()
#    
#    
#    
#    D['bkgSig'] = rhos[parseNumber] 
#    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 16
    
    if parseNumber == 0:
        D['bkgSig'] = 0.005
        D['bkgNo'] = 100
    elif parseNumber == 1:
        D['bkgSig'] == 0.05
        D['bkgNo'] = 100
    elif parseNumber == 2:
        # D['bkgSig'] = .5
        D['bkgNo'] = 100
                
    
    return D