'''
Created on May 26, 2012

@author: dstrauss
'''


import maxwell
import matplotlib.pylab as plt
import scipy.io as io
import numpy as np
import scipy.special as spec

freq = 1e4
P = maxwell.twoDim(freq)
pmc = 31.275936955527357

nx = 99
ny = 99
dx = 5.0
dy = 5.0

eHS = 1.0
sHS = 0.0

P.setspace(nx,ny,dx,dy)
P.setmats(eHS,sHS,ny/2);

P.point_source(49, 49)
# P.makeFD(pmc)
P.setOperators()

P.fwd_solve(0)
M = P.sol[0]

io.savemat('grnP', {'M': M})

epso = 8.854e-12
muo = 4.0*np.pi*1e-7

x = np.array(range(nx))*dx - dx*(nx/2)
print x.shape
print x[0]
print x[98]
Y,X = np.meshgrid(x, x);
dist = np.sqrt(Y**2+X**2) 



c = 1/np.sqrt(muo*epso)

k = (2*np.pi)/(c/freq)
bbx = dx*dx*(1j/4)*spec.hankel1(0,k*dist)
mdpt = nx/2
# the center contains a singularity, but we'll  be forgiving and finite
bbx[mdpt,mdpt] = 0.25*bbx[mdpt+1,mdpt] + 0.25*bbx[mdpt-1,mdpt] + 0.25*bbx[mdpt,mdpt+1] + 0.25*bbx[mdpt,mdpt-1]    

io.savemat('hnkl', {'bbx':bbx})

plt.figure(1)
plt.plot(P.sol[0][:,48:52].real)


plt.figure(2)
plt.subplot(121)
plt.imshow((bbx-P.sol[0]).real)
plt.colorbar()

plt.subplot(122)
plt.imshow((bbx-P.sol[0]).imag)
plt.colorbar()

plt.show()