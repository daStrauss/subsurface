'''
Created on Jun 12, 2012

@author: dstrauss
'''

import forward.flat as flat
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt
import scipy.sparse.linalg as lin

# EE = flat.makeMeA('TE', 1e6, 45*np.pi/180)

EE = flat.makeMeA('TE', 1e6, 45*np.pi/180)    

trm = spio.loadmat('mats/tMat' + repr(1) + '.mat')
pTrue = trm['scrt'].flatten()

EE.initBig(pTrue)

EE.buildROM(100, force=False)


M = (EE.nabla2+EE.getk(1))*EE.Phi

print M.shape
S = np.dot(M.T.conj(),M)
b = np.dot(M.conj().T,EE.rhs)

print 'starting solve'
ur = lin.spsolve(S,b)
print ur.shape



plt.figure(1)
plt.subplot(1,2,1)
plt.imshow(EE.sol[1].reshape(199,199).real)
plt.colorbar()

plt.subplot(1,2,2)
plt.imshow((np.dot(EE.Phi,ur)).reshape(199,199).real)
plt.colorbar()

plt.show()

#MM = flat.makeMeA('TM', 50e3, 45*np.pi/180)
#
#p = np.ones(2000)*0.05
#MM.initBig(p)
#
##nx = 199
##ny = 199
##dx = 5.0
##dy = 5.0
##eHS = 1.0
##sHS = 0.0001
##    
##MM.setspace(nx,ny,dx,dy)
##MM.setOperators()
##MM.makeGradOper()
##MM.setmats(eHS, sHS, ny/2)
##MM.setMs()
##MM.setMd([60,140], [70,95])
##
##p = np.ones(2000)
##MM.sigmap[1] += MM.Md.T*p
#
## rhsz = np.zeros((nx+1,nx+1))
## rhsz[120,120] = -1
#
## rhsx = np.zeros((nx+1,ny))
## rhsy = np.zeros((nx,ny+1))
#
## rhs = np.concatenate((rhsx.flatten(), rhsy.flatten(), rhsz.flatten()))
#
## MM.rhs = rhs
#
#bndl = MM.planeWave()
##MM.fwd_solve(1)
##MM.fwd_solve(0)
#
#fieldsP = MM.parseFields(MM.sol[1])
#fields = MM.parseFields(MM.sol[0])
#
#plt.figure(1)
#plt.imshow((fields[1][20:180,20:180]-bndl[0][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(2)
#plt.imshow((fields[2][20:180,20:180]-bndl[1][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(3)
#plt.imshow((fields[0][20:180,20:180]-bndl[2][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(11)
#plt.imshow((fieldsP[1][20:180,20:180]-bndl[0][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(12)
#plt.imshow((fieldsP[2][20:180,20:180]-bndl[1][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(13)
#plt.imshow((fieldsP[0][20:180,20:180]-bndl[2][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(101)
#plt.imshow((fieldsP[1][20:180,20:180]-fields[1][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(102)
#plt.imshow((fieldsP[2][20:180,20:180]-fields[2][20:180,20:180]).real)
#plt.colorbar()
#
#plt.figure(103)
#plt.imshow((fieldsP[0][20:180,20:180]-fields[0][20:180,20:180]).real)
#plt.colorbar()
#
#
#plt.show()
#
