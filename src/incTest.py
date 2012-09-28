'''
Created on Sep 28, 2012

@author: dstrauss
'''

import forward.flat
import numpy as np
import scipy.io as spio

def test(ica):
    flavor = 'TE'
    freq = 1e6
    incAng = ica*np.pi/180.0
    fwd = forward.flat.makeMeA(flavor, freq, incAng)
    
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    trm = spio.loadmat('mats/tMat' + repr(1) + '.mat')
    pTrue = trm['scrt'].flatten()
    
    fwd.initBig(pTrue, 0.0000)
    
    return fwd


