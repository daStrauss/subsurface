'''
Created on May 25, 2012

@author: dstrauss
'''

import maxwell

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