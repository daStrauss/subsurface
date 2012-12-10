'''
Created on Sep 17, 2012

@author: dstrauss
'''
import numpy as np


D = {'solverType':'sba', 'flavor':'TE', 'numRuns':500, 'expt':'standard', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    bbk = np.int_(np.linspace(10,1600,5))
    romNo, bkg = np.meshgrid(bbk,np.arange(0,100))
    romNo = romNo.flatten()
    bkg = bkg.flatten() 
    
    D['rom'] = romNo[parseNumber]
    D['bkgNo'] = bkg[parseNumber] + 100   
    D['freqs'] = np.array([0.2e3, 6e3, 12e3, 24e3])
    return D