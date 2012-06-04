#!/usr/bin/env python
'''
Created on May 25, 2012

@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import contrastADMM as admm

def bigProj(freq,incAng, ranks):
    S = map(admm.problem, freq, incAng, ranks)
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
    S = map(admm.problem, freq, incAng, ranks)
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

def serial():
    # rho = 1500.0
    rho = 1e-3
    xi = 2e-3
    lmb = 0.0
    uBound = 0.05
    print 'xi = ' + repr(xi) + ' rho = ' + repr(rho)

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
            # run each individual update
            S[ix].runOpt(P)
        
        # aggregate over all
        P = admm.aggregatorSerial(S, lmb, uBound)
        resid[itNo] = np.linalg.norm(P-0.01)

    admm.plotSerial(S, P, resid)


def parallel():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    
#    rho = 1500.0
#    xi = 2e-3
#    lmb = 1e-8
    
    rho = 1e-3
    xi = 2e-3
    lmb = 0.0
    
    
    uBound = 0.05
    print xi
    print rho
    
    #  
    allFreq = np.array([1e3, 3e3, 8e3, 1e4])
    allIncAng = np.ones(allFreq.shape)*45*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = smallProj([allFreq[rank]], [allIncAng[rank]], [rank])
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    S = S[0]
    
    uHat = S.Ms*(S.sol[1].flatten())
    S.initOpt(rho,xi,uHat)
    
    P = np.zeros(S.nRx*S.nRy)
    resid = np.zeros(100)
    
    for itNo in range(2):
        print 'iter no ' + repr(itNo)
        
        S.runOpt(P)
        
        P = admm.aggregatorParallel(S, lmb, uBound, comm)
        resid[itNo] = np.linalg.norm(P-0.01)
        
        
    # do some plotting        
    admm.plotParallel(S,P,resid,rank)
    

    
if __name__ == "__main__":
    parallel()
    
    
    