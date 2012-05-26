'''
Created on May 23, 2012

@author: dstrauss
'''

from maxwell import twoDim
import matplotlib.pyplot as plt
import numpy as np
# if I'm going to be tweaking the class structure, I need to run through
# and redefine what I wa nt here. The re-importing doesn't seem to work.


# from maxwell import findBestAng

# bbx = findBestAng(1e6)

bce = twoDim(1e6)
print 'done got it'
pmc = 22.229964825261945 # best for 1MHz
mc = 1526.4179671752333 # best for 10kHz
nx = 149
ny = nx
dx = 5.0
dy = dx
eHS = 1.0
sHS = 0.00002
sHS = 0.0
bce.setspace(nx,ny,dx,dy)
bce.setmats(eHS,sHS,ny/2);
bce.setMd([55,95],[55,70])
bce.setMs(10)

# bce.point_source(100,50);


# hard code a different sigmap -- interesting
bce.sigmap[1] = (bce.Md.T*np.ones(40*15)*0.001).reshape(nx,ny)

bce.setOperators()

bce.point_source(74, 74)
# bce.te_pw(45*3.141/180)
# bce.point_source(75, 120)
bce.fwd_solve(0)
bce.fwd_solve(1)

bce.setMs(20)
print bce.Ms.size
print bce.Ms.shape

z = bce.Ms*bce.sol[1].flatten('F')

plt.figure(2)
plt.subplot(121)
plt.plot(z.real)
plt.subplot(122)
plt.plot(z.imag)



#
plt.figure(1)
plt.subplot(121)
img = plt.imshow(bce.sol[0].real.T, origin='lower')
# img.set_clim(-1.0,1.0)
plt.title('Real un perturbed')
plt.colorbar()

plt.subplot(122)
img = plt.imshow(bce.sol[1].real.T, origin='lower')
# img.set_clim(-1.0,1.0)
plt.title('Real perturbed')
plt.colorbar()
plt.show()

