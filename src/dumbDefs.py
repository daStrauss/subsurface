'''
Created on Jul 13, 2012

@author: dstrauss
'''

import numpy as np


D = {'solver':'splitField', 'flavor':'TE', 'numRuns':100}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    sigmaBkg, profileNo = np.meshgrid(np.logspace(-3,0,10), np.arange(10))
    sigmaBkg = sigmaBkg.flatten()
    profileNo = profileNo.flatten()
    
    
    D['bkgNo'] = profileNo[parseNumber]
    D['bkgSig'] = sigmaBkg[parseNumber]
    
    return D