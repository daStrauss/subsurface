'''
Created on Jun 12, 2012

@author: dstrauss
'''

import forward.flat as flat
import numpy as np
import scipy.io as spio
import matplotlib.pyplot as plt

# EE = flat.makeMeA('TE', 1e6, 45*np.pi/180)

MM = flat.makeMeA('TM', 1e6, 45*np.pi/180)



nx = 199
ny = 199
dx = 5.0
dy = 5.0
eHS = 1.0
sHS = 0.0
    
MM.setspace(nx,ny,dx,dy)
MM.setOperators()
MM.makeGradOper()
MM.setmats(eHS, sHS, ny/2)
MM.setMs()
MM.setMd([60,140], [70,95])

# rhsz = np.zeros((nx+1,nx+1))
# rhsz[120,120] = -1

# rhsx = np.zeros((nx+1,ny))
# rhsy = np.zeros((nx,ny+1))

# rhs = np.concatenate((rhsx.flatten(), rhsy.flatten(), rhsz.flatten()))

# MM.rhs = rhs

bndl = MM.planeWave()
MM.fwd_solve(0)

ex,ey,hz = MM.parseFields(MM.sol[0])
rex,rey,rhz = MM.parseFields(MM.rhs)


# sex,sey,shz = MM.parseFields(MM.sigmap[0])
