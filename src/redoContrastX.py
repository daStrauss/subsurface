'''
Created on Jul 18, 2012

@author: dstrauss
'''

'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':200, 'expt':'redo'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber    
    return D