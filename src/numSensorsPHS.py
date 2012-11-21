'''
Created on Nov 21, 2012

@author: dstrauss
run the number sensors / number frequencies experiment for the ?phase split alg.
'''
import numpy as np
# import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':4800, 'expt':'noSense', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    noFreqs,noPhis,bkg = np.mgrid[1:7,0:16,0:50]
    noFreqs = noFreqs.flatten()
    noPhis = noPhis.flatten() 
    bkg = bkg.flatten()
    
    allfrq = np.floor(np.linspace(2,150,16))
        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    # D['inc'] = [75.0*np.pi/180.0] -- use defaults!
    D['numSensors'] = allfrq[noPhis[parseNumber]]
    
    D['bkgNo'] = bkg[parseNumber]+100;
    D['incAng'] = np.array([45*np.pi/180.0])
    D['numProcs'] = len(D['freqs'])
    D['rho'] = 1e-3
    D['xi'] = 1e-12
            
            
        
    return D