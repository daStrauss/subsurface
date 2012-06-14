#!/usr/bin/env python
'''
Created on May 25, 2012

@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import scipy.io as spio
import time
import os
import sys

import matplotlib

matplotlib.use('PDF')
# import matplotlib.pyplot as plt

MAXIT = 1

def delegator(solverType, flavor, freq, incAng):
    ''' A function that will allocate the problem instances according to the 'type' given 
    Since I don't mix solvers, it helps to keep the import load low
    '''
    if solverType == 'contrastX':
        import contrastADMM
        S = map(contrastADMM.problem, flavor, freq, incAng)
        return S
    elif solverType == 'splitField':
        import admm
        S = map(admm.problem, freq, incAng, flavor)
        return S
    elif solverType == 'sba':
        import sba
        S = map(sba.problem, flavor, freq, incAng)
        return S
    
def bigProj(S, outDir, testNo):
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    
    trm = spio.loadmat('mats/tMat' + repr(testNo) + '.mat')
  
    pTrue = trm['scrt'].flatten()
    
    for F in S:
        F.fwd.initBig(pTrue)
        F.outDir = outDir
#        
#        F.setspace(nx,ny,dx,dy)
#        F.setmats(eHS,sHS,ny/2);
#        F.setMd([60,140],[70,95])
#        F.setMs(30)   
#        F.setOperators()
#        F.te_pw()
#        F.fwd_solve(0)
#        F.sigmap[1] = F.sigmap[1] + (F.Md.T*pTrue).reshape(nx,ny)
#        
#        F.fwd_solve(1)
#        F.outDir = outDir
    
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
        ti = time.time()
        F.initOpt(uHat,rho,xi,uBound, lmb)
        print 'initalization time ' + repr(time.time()-ti)
    
    # P = np.zeros(80*25)
    P = np.zeros(S[0].nRx*S[0].nRy)
    
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(100):
        
        ti = time.time()
        for ix in range(N):
            # run each individual update
            S[ix].runOpt(P)
        
        print 'iter no ' + repr(itNo) + ' exec time ' + repr(time.time()-ti)
        # aggregate over all
        if solverType == 'sba':
            P += S[0].aggregatorSerial(S)
        else:
            P = S[0].aggregatorSerial(S)
            
        resid[itNo] = np.linalg.norm(P-pTrue)

    S[0].plotSerial(S, P, resid)


def parallel(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1, outDir='basic'):
    '''Parallel solver -- i.e. has MPI calls'''
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    timeFull = time.time()
    
    fout = open('notes'+repr(rank) + '_' +repr(bkgNo) + '.nts', 'w')
    
    fout.write('xi ' + repr(xi) + ' rho = ' + repr(rho) + '\n')
        
    #  
    allFreq = np.array([1e3, 3e3, 13e3, 50e3])
    allIncAng = np.ones(allFreq.shape)*45*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = delegator(solverType, [allFreq[rank]], [allIncAng[rank]], [rank])
    S,pTrue = bigProj(S, outDir, bkgNo)
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    S = S[0]
    
    uHat = S.Ms*(S.sol[1].flatten())
    
    ti = time.time() 
    S.initOpt(uHat, rho, xi, uBound, lmb)
    fout.write('initialization time ' + repr(time.time()-ti) + '\n')
    
    P = np.zeros(S.nRx*S.nRy)
    resid = np.zeros(50)
    
    for itNo in range(50):
        ti = time.time()        
        S.runOpt(P)
        
        # i don't think i can get around this!
        if solverType == 'sba':
            P += S.aggregatorParallel(comm)
        else:
            P = S.aggregatorParallel(comm)
            
        resid[itNo] = np.linalg.norm(P-pTrue)
        fout.write('iter no ' + repr(itNo) + ' exec time = ' + repr(time.time()-ti) + ' rank ' + repr(comm.Get_rank()) +'\n')
        
    # do some plotting        
    S.plotParallel(P,resid,rank)
    S.writeOut()
    
    if rank == 0:
        D = {'Pfinal':P.reshape(S.nRx,S.nRy), 'nProc':nProc, 'resid':resid}
        spio.savemat('pout', D)
        
    fout.write('Solve time = ' + repr(time.time()-timeFull) + '\n')
    fout.close()
        
def semiParallel(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1, outDir='basic'):
    '''semiParallel solver -- i.e. has MPI calls loops locally over different angles of arrival'''
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    timeFull = time.time()
    
    fout = open(outDir + 'notes' + repr(rank) + '_' + repr(bkgNo) + '.nts', 'w')
    
    fout.write('xi ' + repr(xi) + ' rho = ' + repr(rho) + '\n')
        
    #  
    allFreq = np.array([1e3, 3e3, 13e3, 50e3])
#    allIncAng = np.array([75, -75, 45, -45])*np.pi/180
    allIncAng = np.array([75])*np.pi/180
    flavors = ['TE']*allIncAng.size
#     allIncAng = np.ones(allFreq.shape)*45*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = delegator(solverType, flavors, allFreq[rank]*np.ones(allIncAng.shape), allIncAng)
    S,pTrue = bigProj(S, outDir, bkgNo)
    
    N = np.size(S)
    print N
    
    for F in S:
        uHat = F.fwd.Ms*(F.fwd.sol[1].flatten())
        ti = time.time()
        F.initOpt(uHat,rho,xi,uBound, lmb, MAXIT)
        fout.write('initalization time ' + repr(time.time()-ti) + '\n')
    
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    # S = S[0]

    P = np.zeros(S[0].fwd.nRx*S[0].fwd.nRy)
    resid = np.zeros(MAXIT)
    
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
            
        resid[itNo] = np.linalg.norm(P-pTrue)
        fout.write('iter no ' + repr(itNo) + ' exec time = ' + repr(time.time()-ti) + ' rank ' + repr(comm.Get_rank()) +'\n')
        
    # do some plotting        
    for ix in range(N):
        S[ix].plotSemiParallel(P,resid,rank,ix)
        S[ix].writeOut(rank,ix)
    
    if rank == 0:
        D = {'Pfinal':P.reshape(S[0].fwd.nRx,S[0].fwd.nRy), 'nProc':nProc, 'resid':resid}
        spio.savemat(outDir + 'pout_' + solverType + repr(bkgNo), D)
        
    fout.write('Solve time = ' + repr(time.time()-timeFull) + '\n')
    fout.close()
        

    
if __name__ == "__main__":
    # semiParallel('sba', rho=0.005, xi=0.9, uBound=0.05, lmb=0)
    # semiParallel('contrastX')
    semiParallel('splitField', rho=1500, xi =2e-3, uBound = 0.05, lmb = 1e-8, bkgNo = 1)
    # parallel('splitField')
