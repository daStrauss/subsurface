'''
Created on Jun 12, 2012

@author: dstrauss
'''

import te
import tm


def makeMeA(strg, freq, incAng):
    if strg == 'TE':
        return te.solver(freq, incAng, strg)
    elif strg == 'TM':
        return tm.solver(freq, incAng, strg)
    else:
        print 'Easy Tiger: ' + repr(strg) + ' aint around'
    
    