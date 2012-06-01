'''
Created on May 25, 2012

@author: dstrauss
'''

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt

import admm
import scipy.io as io


def bigProj(freq):
    S = map(admm.fieldSplit, freq)
    N = np.size(S)
    nx = 199
    ny = 199
    dx = 5.0
    dy = 5.0
    eHS = 1.0
    sHS = 0.005
    
    
    for F in S:
        F.setspace(nx,ny,dx,dy)
        F.setmats(eHS,sHS,ny/2);
        F.setMd([60,140],[70,95])
        F.setMs(30)   
        F.setOperators()
        F.te_pw(45*np.pi/180)
        F.fwd_solve(0)
        F.sigmap[1] = F.sigmap[1] + (F.Md.T*np.ones(80*25)*0.01).reshape(nx,ny)
        
        F.fwd_solve(1)
    
    return S


def main():
    rho = 1500.0
    xi = 2e-3
    lmb = 1e-8
    uBound = 0.05
    print xi
    print rho

    freq = np.array([1e3, 1e4])
    S = bigProj(freq)
    N = np.size(S)
    
    for F in S:
        uHat = F.Ms*(F.sol[1].flatten())
        
        F.initOpt(rho,xi,F.Ms*(F.sol[1].flatten()))
    
    
    P = np.zeros(80*25)
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(100):
        print 'iter no ' + repr(itNo)
        for ix in range(N):
            S[ix].runOpt(P)
        
        P = admm.aggregateFS(S, lmb, uBound)
        resid[itNo] = np.linalg.norm(P-0.01)
        # io.savemat('iter' + repr(itNo), {'P':P, 'u':S[0].us, 'v':S[0].v, 'E':S[0].E, 'F':S[0].F})
        
        
    plt.figure(383)
    plt.plot(resid)
    
    plt.figure(387)
    plt.imshow(P.reshape(80,25), interpolation='nearest')
    plt.colorbar()
    
    for ix in range(N):
        plt.figure(50+ix)
        vv = S[ix].Ms*S[0].v
        uu = S[ix].Ms*S[0].us
        ub = S[ix].Ms*S[0].ub
        skt = S[ix].uHat-ub
    
        # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})

        plt.plot(np.arange(30), skt.real, np.arange(30), uu.real, np.arange(30), vv.real)
    
    plt.figure(76)
    plt.subplot(121)
    plt.imshow(S[0].us.reshape(199,199).real)
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(S[0].v.reshape(199,199).real)
    plt.colorbar()
    plt.show()
    
    return uHat,ub,uu,S,P

#uHat,ub,us,S,P = solveADMM.main()
    
if __name__ == "__main__":
    main()