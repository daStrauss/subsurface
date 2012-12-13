'''
Created on Dec 11, 2012

@author: dstrauss
'''

import numpy as np
D = {'solverType':'splitField', 'flavor':'TE3D', 'numRuns':125, 'expt':'hugeParm', 'numProcs':1}


def getMyVars(parseNumber, D):
    
    rhos,xis,bkg,frq = np.meshgrid(np.logspace(-4,4,5), np.logspace(-5,-2,5), range(5),np.array([1e3,1e4]) )
    rhos = rhos.flatten()
    bkg = bkg.flatten()
    xis = xis.flatten()
    frq = frq.flatten()
    
    
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  bkg[parseNumber]
    D['numProcs'] = 1
    D['xi'] = xis[parseNumber]
    D['rho'] = rhos[parseNumber]
    D['freqs'] = np.array(frq[parseNumber])
#    D['maxIter'] = 1000
    
    return D 

        