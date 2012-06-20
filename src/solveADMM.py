#!/usr/bin/env python
'''
Created on May 25, 2012

@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import scipy.io as spio
import time
# import os
# import sys

# import matplotlib

# matplotlib.use('PDF')
# import matplotlib.pyplot as plt

MAXIT = 30

def delegator(solverType, flavor, freq, incAng):
    ''' A function that will allocate the problem instances according to the 'type' given 
    Since I don't mix solvers, it helps to keep the import load low
    '''
    if solverType == 'contrastX':
        import contrastADMM
        S = map(contrastADMM.problem, freq, incAng, flavor)
        return S
    elif solverType == 'splitField':
        import admm
        S = map(admm.problem, freq, incAng, flavor)
        return S
    elif solverType == 'sba':
        import sba
        S = map(sba.problem, freq, incAng, flavor)
        return S
    elif solverType == 'biconvex':
        import biconvex
        S = map(biconvex.problem, freq, incAng, flavor)
        return S
    
def bigProj(S, outDir, testNo):
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    
    trm = spio.loadmat('mats/tMat' + repr(testNo) + '.mat')
    pTrue = trm['scrt'].flatten()
    
    for F in S:
        F.fwd.initBig(pTrue)
        F.outDir = outDir
    
    return S,pTrue

def smallProj(S,outDir,testNo):
    '''build for a small project, ie 99x99 '''
    for F in S:
        F.fwd.initSmall(0)
        F.outDir = outDir
    
    pTrue = np.ones((40,10))*0.01
    pTrue = pTrue.flatten()
    return S,pTrue 

def balancingAct(freqs,incAngs,rank,nProc):
    ''' splits the full set of freqs, and incAngs into equal sections according to rank, nProc'''
    allFreqs, allAngs = np.meshgrid(freqs,incAngs)
    allFreqs = allFreqs.flatten()
    allAngs = allAngs.flatten()
    
    nPer = len(allFreqs)/nProc
    assert nPer*nProc == allFreqs.size
    
    lRng = rank*nPer + np.arange(nPer)
    return allFreqs[lRng],allAngs[lRng]
    
    
def semiParallel(solverType, flavor, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1, outDir='basic'):
    '''semiParallel solver -- i.e. has MPI calls loops locally over different angles of arrival'''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    timeFull = time.time()
    
    fout = open(outDir + 'notes' + repr(rank) + '_' + repr(bkgNo) + '.nts', 'w')
    
    fout.write('xi ' + repr(xi) + ' rho = ' + repr(rho) + '\n')
        
    # define the selection that we are interested in
    # remember, they get meshed
    allFreq = np.array([1e3, 3e3, 13e3, 50e3])
    allIncAng = np.array([75, -75, 45, -45])*np.pi/180
    
    # allocate according to the number of processors available
    freqLocal,angLocal = balancingAct(allFreq,allIncAng, rank, nProc)
    
    # switch for local testing
    freqLocal = [freqLocal[2]]; angLocal = [angLocal[2]]
    print freqLocal
    print angLocal
    flavors = [flavor]*len(freqLocal)
    
    # the delegator makes the local set of problems
    S = delegator(solverType, flavors, freqLocal, angLocal)
    
    S,pTrue = smallProj(S, outDir, bkgNo)
    
    N = np.size(S)

    for F in S:
        uHat = F.fwd.Ms*(F.fwd.sol[1].flatten())
        ti = time.time()
        F.initOpt(uHat,rho,xi,uBound, lmb, MAXIT)
        fout.write('initialize time ' + repr(time.time()-ti) + '\n')

    P = np.zeros(S[0].fwd.nRx*S[0].fwd.nRy)
    resid = np.zeros(MAXIT)
    tmvc = np.zeros(MAXIT)
    
    for itNo in range(MAXIT):
        ti = time.time()
        for F in S:        
            objF = F.runOpt(P)
            
            F.obj[itNo] = objF
        
        # i don't think i can get around this!
        if solverType == 'sba':
            P += S[0].aggregatorSemiParallel(S,comm)
        else:
            
            P = S[0].aggregatorSemiParallel(S,comm)
            
        tmvc[itNo] = time.time()-ti    
        resid[itNo] = np.linalg.norm(P-pTrue)
        fout.write('iter no ' + repr(itNo) + ' exec time = ' + repr(time.time()-ti) + ' rank ' + repr(comm.Get_rank()) +'\n')
        fout.flush()
        
    # do some plotting        
    for ix in range(N):
        # S[ix].plotSemiParallel(P,resid,rank,ix)
        S[ix].writeOut(rank,ix)
    
    if rank == 0:
        D = {'Pfinal':P.reshape(S[0].fwd.nRx,S[0].fwd.nRy), 'nProc':nProc, 'resid':resid, \
             'timing':tmvc}
        spio.savemat(outDir + 'pout_' + solverType + repr(bkgNo), D)
        
    fout.write('Solve time = ' + repr(time.time()-timeFull) + '\n')
    fout.close()
        

    
if __name__ == "__main__":
#    semiParallel('sba', 'TE', rho=0.005, xi=0.9, uBound=0.05, lmb=0)
    # semiParallel('biconvex', 'TE', rho=0.001, xi=1e-5, lmb=0, uBound=0.05,bkgNo=1)
    semiParallel('contrastX', 'TM', rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1)
#    semiParallel('splitField','TE', rho=1500, xi =2e-3, uBound = 0.05, lmb = 1e-8, bkgNo = 1)
    # parallel('splitField')
    
