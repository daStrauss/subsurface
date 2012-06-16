'''
Created on Jun 12, 2012

@author: dstrauss
'''

import te
import tm


def makeMeA(strg, freq, incAng):
    if strg == 'TE':
        return te.solver(freq, incAng)
    elif strg == 'TM':
        return tm.solver(freq, incAng)
    else:
        print 'Easy Tiger: ' + repr(strg) + ' aint around'
    
    