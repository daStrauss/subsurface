'''
Created on Jun 14, 2012

@author: dstrauss
'''

import forward.flat

class optimizer(object):
    '''base class for optimization routines'''
    
    def __init__(self,freq,incAng,flavor):
        self.fwd = forward.flat.makeMeA(flavor, freq, incAng)
        