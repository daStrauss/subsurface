#!/usr/bin/env python
'''
Created on May 25, 2012
Copyright © 2013
The Board of Trustees of The Leland Stanford Junior University.
All Rights Reserved

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.


@author: dstrauss
'''

from mpi4py import MPI
import numpy as np
import scipy.io as spio
import time
import solverDefaults
import importlib
# import os
# import sys

# import matplotlib

# matplotlib.use('PDF')
# import matplotlib.pyplot as plt


def delegator(solverType, flavor, freq, incAng):
    ''' A function that will allocate the problem instances according to the 'type' given 
    Since I don't mix solvers, it helps to keep the import load low
    '''
    
    ''' change to use a more direct importing method '''
    
    optm = importlib.import_module('slvr.' + solverType)
    S = map(optm.problem, freq,incAng,flavor)
    
    return S
    
    
    
def bigProj(S, D):
    ''' Define a big project, with a tag and a test No -- will draw from ../mats'''
    trm = spio.loadmat('mats/tMat' + repr(D['bkgNo']+1) + '.mat')
    pTrue = trm['scrt'].flatten()
    
    for F in S:
        F.fwd.initBig(pTrue, bkgSig=D['bkgSig'], numSensors=D['numSensors'])
        F.outDir = D['outDir']
        if not D['rom']:
            F.fwd.isRom = False
        else:
            F.fwd.buildROM(D['rom'],force=True)
            
    return S,pTrue

def smallProj(S,D):
    '''build for a small project, ie 99x99 '''
    pTrue = np.ones((40,10))*0.01
    pTrue = pTrue.flatten()
    
    for F in S:
        F.fwd.initSmall(pTrue,D['bkgSig'])
        F.outDir = D['outDir']
        if not D['rom']:
            F.fwd.isRom = False
        else:
            F.fwd.buildROM(D['rom'],force=True)
    
    return S,pTrue 

def feas3dProj(S,D):
    '''build for a small project, ie 41x41x41 '''
    # pTrue = np.ones((14,7,14))*0.01
    # pTrue = np.ones((7,3,7))*0.005
    # pTrue[:3,:,:] = 0.02
    # pTrue = pTrue.flatten()
    trm = spio.loadmat('mats/t3dMat' + repr(D['bkgNo']+1) + '.mat')
    pTrue = trm['scrt'].flatten()
    
    for F in S:
        F.fwd.init3Dfeas(pTrue,D['bkgSig'])
        F.outDir = D['outDir']
        if not D['rom']:
            F.fwd.isRom = False
        else:
            F.fwd.buildROM(D['rom'],force=True)
    
    return S,pTrue 

def balancingAct(D,rank,nProc):
    ''' splits the full set of freqs, and incAngs into equal sections according to rank, nProc'''
    
    allFreqs, allAngs = np.meshgrid(D['freqs'], D['inc'])
    allFreqs = allFreqs.flatten()
    allAngs = allAngs.flatten()
    # print(D['flavor'])
    
    if D['flavor'] == 'both':
        allFlav = ['TE']*len(allFreqs) + ['TM']*len(allFreqs)
        allFreqs = np.tile(allFreqs,2)
        allAngs = np.tile(allAngs,2)
    else:
        allFlav = [D['flavor']]*len(allFreqs)
    
    # print(repr(allFreqs.shape))
    # print(repr(nProc))
    
    nPer = len(allFreqs)/nProc
    assert nPer*nProc == allFreqs.size
    
    # print(allFlav)
    lRng = rank*nPer + np.arange(nPer)
    # print(allFlav[lRng] + ' ' + repr(lRng))
    return allFreqs[lRng],allAngs[lRng],[allFlav[lRng]]
    
    
def semiParallel(solverType, flavor, **kwargs): 
    ''' The main coordination routine'''
    D = solverDefaults.getDefaults(solverType, flavor, kwargs)
    
    '''semiParallel solver -- i.e. has MPI calls loops locally over different angles of arrival'''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    timeFull = time.time()
    
    fout = open(D['outDir'] + 'notes' + repr(rank) + '_' + repr(D['ix']) + '.nts', 'w')
    
    fout.write('xi ' + repr(D['xi']) + ' rho = ' + repr(D['rho']) + '\n')

    
    # allocate according to the number of processors available
    freqLocal,angLocal,flavLocal = balancingAct(D, rank, nProc)
    
    fout.write('rank ' + repr(rank) + ' ' + repr(flavLocal) + ' frq ' + repr(freqLocal))
    # switch for local testing
    # freqLocal = [freqLocal[2]]; angLocal = [angLocal[2]]
    
    # the delegator makes the local set of problems
    S = delegator(solverType, flavLocal, freqLocal, angLocal)
    print flavLocal
    
    if flavLocal[0] == 'TE3D':
        print 'is this happening'
        S,pTrue = feas3dProj(S,D)
    else:
        print 'or am I not making it'
        S,pTrue = bigProj(S, D)
    
    N = np.size(S)

    for F in S:
        uHat = F.fwd.Ms*(F.fwd.sol[1].flatten())
        
        noi = np.random.randn(F.fwd.nSen) + 1j*np.random.randn(F.fwd.nSen)
        noi = D['relNoi']*(np.linalg.norm(uHat)/np.linalg.norm(noi))*noi
        
        uHat = uHat + noi
        ti = time.time()
        F.initOpt(uHat,D)
        fout.write('initialize time ' + repr(time.time()-ti) + '\n')
    
    
    P = np.zeros(S[0].fwd.getPSize())
    resid = np.zeros(D['maxIter'])
    tmvc = np.zeros(D['maxIter'])
    

    for itNo in range(D['maxIter']):
        ti = time.time()
        for F in S:        
            objF = F.runOpt(P)
            
            F.obj[itNo] = objF
        print 'rank ' + repr(rank) + ' ' + repr(flavLocal) +'finished independent ' + repr(itNo)
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
        D['Pfinal'] = S[0].fwd.parseP(P)
        D['nProc'] = nProc
        D['resid'] = resid
        D['timing'] = tmvc
        D['pTrue'] = pTrue
        spio.savemat(D['outDir'] + 'pout_' + solverType + repr(D['ix']), D)

        
    fout.write('Solve time = ' + repr(time.time()-timeFull) + '\n')
    fout.close()
        

    
if __name__ == "__main__":
#    semiParallel('sba', 'TM', rho=0.005, xi=0.9, uBound=0.05, lmb=0,bkgNo=1)
    # semiParallel('biconvex', 'TM')
#    semiParallel('sba', 'TE', freqs=[1e3], inc=[45*np.pi/180], maxIter=100, rom=75, lmb=0)
    semiParallel('sba', 'TE', maxIter=1000, rom=75, lmb=0)
#    semiParallel('contrastX', 'TM', rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1)
#    semiParallel('splitField','TE', rho=1500, xi =2e-3, uBound = 0.05, lmb = 1e-8, bkgNo = 1)

