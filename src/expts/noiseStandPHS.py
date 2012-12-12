'''
Created on Dec 10, 2012

@author: dstrauss
'''
import numpy as np

D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':1000, 'expt':'bringDa', 'numProcs':4}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    
    snr,bkg = np.meshgrid(np.logspace(-5,0,20),range(50))
        
    snr = snr.flatten()
    bkg = bkg.flatten()                 

    D['relNoi'] = snr[parseNumber]        
    D['bkgNo'] = bkg[parseNumber]+100
    
    D['inc'] = np.array([75.0*np.pi/180.0])
    D['numProcs'] = 4
    D['rho'] = 1e-3
    D['xi'] = 1e-12

    
        
    return D