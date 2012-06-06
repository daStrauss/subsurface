#!/usr/bin/env python
'''
Created on May 25, 2012

@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import scipy.io as spio


def delegator(solverType, freq, incAng, ranks):
    ''' A function that will allocate the problem instances according to the 'type' given '''
    if solverType == 'contrastX':
        import contrastADMM
        S = map(contrastADMM.problem, freq, incAng, ranks)
        return S
    elif solverType == 'splitField':
        import admm
        S = map(admm.problem, freq, incAng, ranks)
        return S
    elif solverType == 'sba':
        import sba
        S = map(sba.problem, freq, incAng, ranks)
        return S
    
def bigProj(S,tag, testNo):
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    
    nx = 199
    ny = 199
    dx = 5.0
    dy = 5.0
    eHS = 1.0
    sHS = 0.005
    
    F = spio.loadmat('../mats/tMat' + repr(testNo) + '.mat')
    pTrue = F['scrt'].flatten()
    
    for F in S:
        F.setspace(nx,ny,dx,dy)
        F.setmats(eHS,sHS,ny/2);
        F.setMd([60,140],[70,95])
        F.setMs(30)   
        F.setOperators()
        F.te_pw()
        F.fwd_solve(0)
        F.sigmap[1] = F.sigmap[1] + (F.Md.T*pTrue).reshape(nx,ny)
        
        F.fwd_solve(1)
        
        F.tag = '_' + repr(testNo) + tag
    
    return S,pTrue

def smallProj(S):
    '''build for a small project, ie 99x99 '''
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

def serial(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1):

    # rho = 1500.0
    # ho = 1e-3
    # xi = 2e-3
    # lmb = 0.0
    
    print 'xi = ' + repr(xi) + ' rho = ' + repr(rho)

    freq = np.array([1e3, 1e4])
    incAng = np.array([45.0, 45.0])*np.pi/180.0
    ranks = np.arange(np.size(freq))
    
    S = delegator(solverType, freq, incAng, ranks)
    # S = smallProj(S)
    S,pTrue = bigProj(S,'basic',bkgNo)
    N = np.size(S)
    
    for F in S:
        uHat = F.Ms*(F.sol[1].flatten())
        
        F.initOpt(uHat,rho,xi,uBound, lmb)
    
    
    # P = np.zeros(80*25)
    P = np.zeros(S[0].nRx*S[0].nRy)
    
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(100):
        print 'iter no ' + repr(itNo)
        for ix in range(N):
            # run each individual update
            S[ix].runOpt(P)
        
        # aggregate over all
        if solverType == 'sba':
            P += S[0].aggregatorSerial(S)
        else:
            P = S[0].aggregatorSerial(S)
            
        resid[itNo] = np.linalg.norm(P-pTrue)

    S[0].plotSerial(S, P, resid)


def parallel(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1):
    '''Parallel solver -- i.e. has MPI calls'''
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    
    uBound = 0.05
    print xi
    print rho
    
    #  
    allFreq = np.array([1e3, 3e3, 8e3, 1e4])
    allIncAng = np.ones(allFreq.shape)*45*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = delegator(solverType, [allFreq[rank]], [allIncAng[rank]], [rank])
    S,pTrue = bigProj(S, 'basic', bkgNo)
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    S = S[0]
    
    uHat = S.Ms*(S.sol[1].flatten())
    
    S.initOpt(uHat, rho, xi, uBound, lmb)
    
    P = np.zeros(S.nRx*S.nRy)
    resid = np.zeros(1000)
    
    for itNo in range(1000):
        print 'iter no ' + repr(itNo)
        
        S.runOpt(P)
        
        # i don't think i can get around this!
        if solverType == 'sba':
            P += S.aggregatorParallel(comm)
        else:
            P = S.aggregatorParallel(comm)
            
        resid[itNo] = np.linalg.norm(P-pTrue)
        
        
    # do some plotting        
    S.plotParallel(P,resid,rank)
    S.writeOut()
    
    if rank == 0:
        D = {'Pfinal':P.reshape(S.nRx,S.nRy), 'nProc':nProc, 'resid':resid}
        spio.savemat('pout', D)
        
        

    
if __name__ == "__main__":
    # parallel('sba', rho=0.005, xi=0.9, uBound=0.05, lmb=0)
    # parallel('contrastX')
    # parallel('splitField', rho=1500, xi =2e-4, uBound = 0.05, lmb = 1e-8, bkgNo = 1)
    serial('splitField')
