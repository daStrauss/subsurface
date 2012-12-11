'''
Created on Jul 14, 2012

@author: dstrauss
'''

'''
Created on Jul 13, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'sba', 'flavor':'TE', 'numRuns':1100, 'expt':'standard', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber    
    return D