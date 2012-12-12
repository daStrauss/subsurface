'''
Created on Dec 11, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'splitField', 'flavor':'TE3D', 'numRuns':125, 'expt':'hugeParm', 'numProcs':1}


def getMyVars(parseNumber, D):
    
    rhos,xis,frq = np.meshgrid(np.logspace(-4,4,5), np.logspace(-10,-4,5), np.logspace(3,5,5))
    rhos = rhos.flatten()
    frq = frq.flatten()
    xis = xis.flatten()
    
    
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['numProcs'] = 1
    D['xi'] = xis[parseNumber]
    D['rho'] = rhos[parseNumber]
    D['freqs'] = np.array(frq[parseNumber])
#    D['maxIter'] = 1000
    
    return D 

        