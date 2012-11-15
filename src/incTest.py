'''
Created on Sep 28, 2012

@author: dstrauss
'''

import forward.flat
import numpy as np
import scipy.io as spio
import time
from superSolve import wrapCvxopt

def test(ica):
    flavor = 'TE3D'
    freq = 1e4
    incAng = ica*np.pi/180.0
    
    strt = time.time()
    fwd = forward.flat.makeMeA(flavor, freq, incAng)
    
    bkgSig=0.005
    numSensors=4
    
    fwd.setspace(21,21,21,5.0,5.0,5.0)
    fwd.setmats(1,bkgSig,21/2)
    fwd.setOperators()
    fwd.makeGradOper()
    
    fwd.setMs(numSensors)
    # for 31 use:
    # fwd.setMd([8,22], [8,15], [8,22])
    # for 21 use:
    # fwd.setMd([7,13], [7,9], [7,13])
    
    
    # fwd.sigmap[1] += self.Md.T*p
    # fwd.pointSource(10,10,10)
#        self.planeWave()
    fwd.planeWave()
    fwd.fwd_solve(0)
    
    A = fwd.nabla2+fwd.getk(0)
    sol = wrapCvxopt.linsolve(A, fwd.rhs)
    
    # fwd.fwd_solve(0)
    # self.fwd_solve(1)
    print 'Run time = ' + repr(time.time()-strt)
    [ex,ey,ez] = fwd.parseFields(sol)
    rhx,rhy,rhz = fwd.parseFields(fwd.rhs)
    
    
    D = {'ox':fwd.nabla2+fwd.getk(0),'rhs':fwd.rhs,'nbl':fwd.nabla2,'ex':ex, 'ey':ey, 'ez':ez,\
         'rhx':rhx,'rhy':rhy,'rhz':rhz}
    
    spio.savemat('threeDTest', D)
    
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    # trm = spio.loadmat('mats/tMat' + repr(1) + '.mat')
    # pTrue = trm['scrt'].flatten()
    
    # fwd.initBig(pTrue, 0.0000)
    
    return fwd

if __name__=='__main__':
    test(45)
