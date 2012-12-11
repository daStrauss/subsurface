'''
Created on Aug 30, 2012

@author: dstrauss
'''


import numpy as np
# import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':3600, 'expt':'incConds', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    if parseNumber < 1800:
        noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,0:50]
        noFreqs = noFreqs.flatten()
        noPhis = noPhis.flatten() 
        bkg = bkg.flatten()
        lcp = parseNumber
    else:
        noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,50:100]
        noFreqs = noFreqs.flatten()
        noPhis = noPhis.flatten() 
        bkg = bkg.flatten()
        lcp = parseNumber-1800
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[lcp]))
    D['inc'] = (np.linspace(-75,75,noPhis[lcp])*np.pi/180.0)
    D['bkgNo'] = bkg[lcp]+100
    D['numProcs'] = len(D['freqs'])*len(D['inc'])
    
        
    return D