'''
Created on Sep 17, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':100, 'expt':'fullFreqs', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber + 100   
    D['freqs'] = np.logspace(2.3010, 4.6990,16) 
    D['inc'] = np.array([75.0])*np.pi/180
    return D