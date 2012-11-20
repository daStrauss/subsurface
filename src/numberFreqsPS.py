'''
Created on Nov 20, 2012

@author: dstrauss

extend the number frequencies simulations to the phase split method
ok so this is awkward. the ability to incorporate multiple frequencies is limited.

'''


import numpy as np
# import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]

# F = {'goTo':np.arange(1800)}

D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':1800, 'expt':'incConds', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,0:50]
    noFreqs = noFreqs.flatten()
    noPhis = noPhis.flatten() 
    bkg = bkg.flatten()
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    D['inc'] = (np.linspace(-75,75,noPhis[parseNumber])*np.pi/180.0)
    D['bkgNo'] = bkg[parseNumber]+100
    D['numProcs'] = len(D['freqs'])*len(D['inc'])
    
        
        
    return D