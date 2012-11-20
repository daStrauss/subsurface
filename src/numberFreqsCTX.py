'''
Created on Aug 30, 2012

@author: dstrauss
'''


import numpy as np
# import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]

F = {'goTo':np.arange(1800)}

D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':3600, 'expt':'incConds', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    if parseNumber < 1800:
        noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,0:50]
        noFreqs = noFreqs.flatten()
        noPhis = noPhis.flatten() 
        bkg = bkg.flatten()
    else:
        noFreqs,noPhis,bkg = np.mgrid[1:7,1:7,50:100]
        noFreqs = noFreqs.flatten()
        noPhis = noPhis.flatten() 
        bkg = bkg.flatten()
        
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    D['inc'] = (np.linspace(-75,75,noPhis[parseNumber])*np.pi/180.0)
    D['bkgNo'] = bkg[parseNumber]+100
    D['numProcs'] = len(D['freqs'])*len(D['inc'])
    
    if parseNumber in F['goTo']:
        print 'here we go'
    else:
        D['numProcs'] = 0
        
        
    return D