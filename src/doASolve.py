'''
Created on May 23, 2012

@author: dstrauss
'''

from maxwell import twoDim
# if I'm going to be tweaking the class structure, I need to run through
# and redefine what I wa nt here. The re-importing doesn't seem to work.


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
bce.sigmap[1][55:95,55:70] = 0.001
print 'sigmap'
bce.make_operators(pmc)
print 'mkop'
# bce.te_pw(45*3.141/180)
bce.point_source(75, 120)
print 'pw'
bce.fwd_solve(0)
bce.fwd_solve(1)


        
    
    
