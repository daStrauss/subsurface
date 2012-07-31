'''
Created on Jul 23, 2012

@author: dstrauss
'''


import numpy as np
import scipy.io as spio

F = spio.loadmat('incCondGo.mat')
numRuns = F['goTo'].shape[0]

D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':3600, 'expt':'incConds', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,0:100]
    noFreqs = noFreqs.flatten()
    noPhis = noPhis.flatten() 
    bkg = bkg.flatten()
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    D['inc'] = np.round(np.linspace(-75,75,noPhis[parseNumber]))
    D['bkgNo'] = bkg[parseNumber]
    D['numProcs'] = len(D['freqs'])*len(D['inc'])
    
    if parseNumber in F['goTo']:
        print 'here we go'
    else:
        D['numProcs'] = 0
        
        
    return D