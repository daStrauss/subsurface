'''
Created on Nov 29, 2012

@author: dstrauss
'''
import numpy as np
D = {'solverType':'ctxdev', 'flavor':'TE3D', 'numRuns':125, 'expt':'hugeParm', 'numProcs':1}


def getMyVars(parseNumber, D):
    
    rhos,reg,frq = np.meshgrid(np.logspace(-4,4,5), np.logspace(-10,-4,5), np.logspace(3,5,5))
    rhos = rhos.flatten()
    frq = frq.flatten()
    reg = reg.flatten()
    
    
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['numProcs'] = 1
    D['reg'] = reg[parseNumber]
    D['rho'] = rhos[parseNumber]
    D['freqs'] = np.array(frq[parseNumber])
    
    return D 

        
