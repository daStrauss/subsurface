'''
Created on May 25, 2012

@author: dstrauss
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import admm



def bigProj(freq):
    S = map(admm.fieldSplit, freq)
    N = np.size(S)
    nx = 149
    ny = 149
    dx = 5.0
    dy = 5.0
    eHS = 1.0
    sHS = 0.005
    
    
    for ix in range(N):
        S[ix].setspace(nx,ny,dx,dy)
        S[ix].setmats(eHS,sHS,ny/2);
        S[ix].setMd([55,95],[55,70])
        S[ix].setMs(10)   
        S[ix].setOperators()
        S[ix].te_pw(45*3.141/180)
        S[ix].fwd_solve(0)
        S[ix].sigmap[1] = (S[ix].Md.T*np.ones(40*15)*0.01).reshape(nx,ny)
        
        S[ix].fwd_solve(1)
    
    return S


def main():
    rho = 1500.0
    xi = 2e-3
    lmb = 0.0
    uBound = 0.05

    freq = np.array([1e4])
    S = bigProj(freq)
    N = np.size(S)
    
    for ix in range(N):
        uHat = S[ix].Ms*S[ix].sol[1].flatten()
        S[ix].initOpt(rho,xi,uHat)
        
    P = np.zeros(40*15)
    resid = np.zeros(100)
    
    for itNo in range(100):
        print 'iter no ' + repr(itNo)
        for ix in range(N):
            S[ix].runOpt(P)
        
        P = admm.aggregateFS(S, lmb, uBound)
        resid[itNo] = np.linalg.norm(P-0.01)
        
    plt.figure(383)
    plt.plot(resid)
    
    plt.figure(387)
    plt.imshow(P.reshape(40,15))
    plt.colorbar()
    plt.show()
    
    
if __name__ == "__main__":
    main()