'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':1100, 'expt':'standard'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber    
    return D