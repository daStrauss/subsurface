'''
Created on May 25, 2012

@author: dstrauss
'''
import numpy as np
import admm
import solveADMM

rho = 1500.0
xi = 2e-3
lmb = 0.0

freq = np.array([1e4])
S = solveADMM.bigProj(freq)
N = np.size(S)

for ix in range(N):
    uHat = S[ix].Ms*S[ix].sol[1].flatten()
    S[ix].initOpt(rho,xi,uHat)
    
P = np.zeros(40*15)
resid = np.zeros(100)