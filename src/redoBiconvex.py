'''
Created on Jul 19, 2012

@author: dstrauss
'''

D = {'solverType':'biconvex', 'flavor':'TE', 'numRuns':200, 'expt':'standard'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber    
    return D