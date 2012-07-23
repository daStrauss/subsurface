'''
Created on Jul 23, 2012

@author: dstrauss
'''


import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':3600, 'expt':'standard', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    noFreqs = noFreqs.flatten()
    noPhis = noPhis.flatten() 
    bkg = bkg.flatten()
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    D['inc'] = np.round(np.linspace(-180,180,noPhis[parseNumber]))
    D['bkgNo'] = bkg[parseNumber]
    D['numProcs'] = len(D['freqs'])*len(D['inc'])
    
        
    return D