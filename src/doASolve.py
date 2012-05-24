'''
Created on May 23, 2012

@author: dstrauss
'''

from maxwell import twoDim
import matplotlib.pyplot as plt
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
bce.setspace(nx,ny,dx,dy)
print 'setspace'
# bce.point_source(100,50);
bce.setmats(eHS,sHS,ny/2);
print 'setmats'

# hard code a different sigmap -- interesting
bce.sigmap[1][55:95,55:70] = 0.001
print 'sigmap'
bce.setOperators()
print 'mkop'
bce.te_pw(45*3.141/180)
# bce.point_source(75, 120)
print 'pw'
bce.fwd_solve(0)
bce.fwd_solve(1)


plt.figure(1)
plt.subplot(121)
plt.imshow(bce.sol[0].real)
plt.title('Real part')
plt.colorbar()

plt.subplot(122)
plt.imshow(bce.sol[0].imag)
plt.title('Imag part')
plt.colorbar()
plt.show()


