'''
Created on Sep 28, 2012

@author: dstrauss
'''

import forward.flat
import numpy as np
import scipy.io as spio

def test(ica):
    flavor = 'TE3D'
    freq = 1e7
    incAng = ica*np.pi/180.0
    
    
    fwd = forward.flat.makeMeA(flavor, freq, incAng)
    
    bkgSig=0.0
    numSensors=10
    
    fwd.setspace(25,25,25,5.0,5.0,5.0)
    fwd.setmats(1,bkgSig,21/2)
    fwd.setOperators()
    fwd.makeGradOper()
    # fwd.setMs(numSensors)
    # fwd.setMd([60,140], [70,95])
    # fwd.sigmap[1] += self.Md.T*p
    fwd.pointSource(10,12,10)
#        self.planeWave()
    fwd.fwd_solve(0)
    # self.fwd_solve(1)
    [ex,ey,ez] = fwd.parseFields(fwd.sol[0])
    
    D = {'ex':ex, 'ey':ey, 'ez':ez}
    
    spio.savemat('threeDTest', D)
    
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    # trm = spio.loadmat('mats/tMat' + repr(1) + '.mat')
    # pTrue = trm['scrt'].flatten()
    
    # fwd.initBig(pTrue, 0.0000)
    
    return fwd


