'''
Created on Dec 10, 2012

@author: dstrauss
'''
import numpy as np

D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':2000, 'expt':'bringDa', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    
    snr,bkg = np.meshgrid(np.logspace(-5,0,20),range(100))
        
    snr = snr.flatten()
    bkg = bkg.flatten()                 

    D['relNoi'] = snr[parseNumber]        
    D['bkgNo'] = bkg[parseNumber]+100
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), 16))
    
    D['inc'] = np.array([75.0*np.pi/180.0])
    D['numProcs'] = 16
    
    D['rho'] = 7.7e-4
    D['xi'] = 6.0e-12
    
#    D['rho'] = 1e-3
#    D['xi'] = 1e-12

    
        
    return D