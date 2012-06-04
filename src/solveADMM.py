#!/usr/bin/env python
'''
Created on May 25, 2012

@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt

import admm

def bigProj(freq,incAng, ranks):
    S = map(admm.fieldSplit, freq, incAng, ranks)
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
        F.te_pw()
        F.fwd_solve(0)
        F.sigmap[1] = F.sigmap[1] + (F.Md.T*np.ones(80*25)*0.01).reshape(nx,ny)
        
        F.fwd_solve(1)
    
    return S

def smallProj(freq, incAng, ranks):
    '''build for a small project, ie 99x99 '''
    S = map(admm.fieldSplit, freq, incAng, ranks)
    nx = 99
    ny = 99
    dx = 5.0
    dy = 5.0
    eHS = 1.0
    sHS = 0.005
    
    for F in S:
        F.setspace(nx,ny,dx,dy)
        F.setmats(eHS,sHS,ny/2)
        F.setMd([30, 70], [35, 45])
        F.setMs(15)
        F.setOperators()
        F.te_pw()
        F.fwd_solve(0)
        F.sigmap[1] = F.sigmap[1] +(F.Md.T*np.ones(40*10)*0.01).reshape(nx,ny)
        
        F.fwd_solve(1)
    return S

def single():
    rho = 1500.0
    xi = 2e-3
    lmb = 1e-8
    uBound = 0.05
    print xi
    print rho

    freq = np.array([1e3, 1e4])
    incAng = np.array([45.0, 45.0])*np.pi/180.0
    ranks = np.arange(np.size(freq))
    
    S = smallProj(freq, incAng, ranks)
    N = np.size(S)
    
    for F in S:
        uHat = F.Ms*(F.sol[1].flatten())
        
        F.initOpt(rho,xi,F.Ms*(F.sol[1].flatten()))
    
    
    # P = np.zeros(80*25)
    P = np.zeros(S[0].nRx*S[0].nRy)
    
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(2):
        print 'iter no ' + repr(itNo)
        for ix in range(N):
            S[ix].runOpt(P)
        
        P = admm.aggregateFS(S, lmb, uBound)
        resid[itNo] = np.linalg.norm(P-0.01)
        # io.savemat('iter' + repr(itNo), {'P':P, 'u':S[0].us, 'v':S[0].v, 'E':S[0].E, 'F':S[0].F})
        
        
    plt.figure(383)
    plt.plot(resid)
    
    plt.figure(387)
    plt.imshow(P.reshape(S[0].nRx,S[0].nRy), interpolation='nearest')
    plt.colorbar()
    
    for ix in range(N):
        plt.figure(50+ix)
        vv = S[ix].Ms*S[0].v
        uu = S[ix].Ms*S[0].us
        ub = S[ix].Ms*S[0].ub
        skt = S[ix].uHat-ub
    
        # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})

        plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real, np.arange(S[0].nSen), vv.real)
    
    plt.figure(76)
    plt.subplot(121)
    plt.imshow(S[0].us.reshape(S[0].nx,S[0].ny).real)
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(S[0].v.reshape(S[0].nx,S[0].ny).real)
    plt.colorbar()
    plt.show()
    
#     S[0].writeOut(1,P)
    
    return uHat,ub,uu,S,P

#uHat,ub,us,S,P = solveADMM.main()


def parallel():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
    rho = 1500.0
    xi = 2e-3
    lmb = 1e-8
    uBound = 0.05
    print xi
    print rho

    allFreq = np.array([1e3, 1e4])
    allIncAng = np.array([45.0, 45.0])*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = smallProj([allFreq[rank]], [allIncAng[rank]], [rank])
    N = np.size(S)
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    S = S[0]
    
    uHat = S.Ms*(S.sol[1].flatten())
    S.initOpt(rho,xi,uHat)
    
    # P = np.zeros(80*25)
    P = np.zeros(S.nRx*S.nRy)
    
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(2):
        print 'iter no ' + repr(itNo)
        
        S.runOpt(P)
        
        P = admm.aggregateParallel(S, lmb, uBound, comm)
        resid[itNo] = np.linalg.norm(P-0.01)
        # io.savemat('iter' + repr(itNo), {'P':P, 'u':S[0].us, 'v':S[0].v, 'E':S[0].E, 'F':S[0].F})
        
    if rank==0:
        # then print some figures   
        plt.figure(383)
        plt.plot(resid)
    
        plt.figure(387)
        plt.imshow(P.reshape(S.nRx,S.nRy), interpolation='nearest')
        plt.colorbar()
    
#        for ix in range(N):
#            plt.figure(50+ix)
        vv = S.Ms*S.v
        uu = S.Ms*S.us
        ub = S.Ms*S.ub
        skt = S.uHat-ub
    
        # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})
        plt.figure(40)
        plt.plot(np.arange(S.nSen), skt.real, np.arange(S.nSen), uu.real, np.arange(S.nSen), vv.real)
    
        plt.figure(76)
        plt.subplot(121)
        plt.imshow(S.us.reshape(S.nx,S.ny).real)
        plt.colorbar()
    
        plt.subplot(122)
        plt.imshow(S.v.reshape(S.nx,S.ny).real)
        plt.colorbar()
        plt.show()
    

    
if __name__ == "__main__":
    parallel()
    
    
    