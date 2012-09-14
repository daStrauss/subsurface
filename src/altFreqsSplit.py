'''
Created on Sep 14, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':100, 'expt':'altFreq'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber + 100   
    D['freqs'] = np.array([0.2e3, 6e3, 12e3, 24e3])
    return D