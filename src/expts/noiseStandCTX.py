'''
Created on Dec 10, 2012

@author: dstrauss
'''

'''
Created on Aug 30, 2012

@author: dstrauss
'''


import numpy as np

D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':2000, 'expt':'bringDa', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    
    snr,bkg = np.meshgrid(np.logspace(-5,0,20),range(1000))
        
    snr = snr.flatten()
    bkg = bkg.flatten()                 

    D['relNoi'] = snr[parseNumber]        
    D['bkgNo'] = bkg[parseNumber]+100

    
        
    return D