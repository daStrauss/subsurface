'''
Created on Jun 12, 2012

@author: dstrauss
'''

import forward.flat as flat
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt

EE = flat.makeMeA('TE', 1e4, 45*np.pi/180)

MM = flat.makeMeA('TM', 1e4, 45*np.pi/180)



nx = 199
ny = 199
dx = 5.0
dy = 5.0
eHS = 1.0
sHS = 0.005
    
MM.setspace(nx,ny,dx,dy)
MM.setOperators()
MM.makeGradOper()
MM.setmats(eHS, sHS, ny/2)
MM.setMs()
MM.setMd([60,140], [70,95])

print MM.Ms.shape
    
outDir = './'
testNo = 1
F = spio.loadmat('mats/tMat' + repr(testNo) + '.mat')
  
pTrue = F['scrt'].flatten()
    
EE.setspace(nx,ny,dx,dy)
EE.setmats(eHS,sHS,ny/2);
EE.setMd([60,140],[70,95])
EE.setMs(30)   
EE.setOperators()
EE.makeGradOper()
EE.planeWave()
EE.fwd_solve(0)
EE.sigmap[1] = EE.sigmap[1] + (EE.Md.T*pTrue).reshape(nx,ny)
        
EE.fwd_solve(1)
EE.outDir = outDir

#plt.figure(1)
#plt.imshow(EE.sol[0].real)
#plt.colorbar()
#
#plt.figure(2)
#plt.imshow(EE.sol[1].real)
#plt.colorbar()
#plt.title('with scatter')
#
#plt.show()



