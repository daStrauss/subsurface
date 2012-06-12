'''
Created on Jun 12, 2012

@author: dstrauss
'''

import te
import tm


def makeMeA(strg):
    if strg == 'TE':
        return te.solver(5,99,99)
    elif strg == 'TM':
        return tm.solver(5,99,99)
    
    