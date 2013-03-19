'''
Created on Mar 19, 2013

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TM', 'numRuns':1100, 'expt':'standard', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber    
    return D