'''
Created on Oct 22, 2012

@author: dstrauss

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

Implementation of the projection operator method for solving the subsurface imaging problem
'''

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import sparseTools as spt
import scipy.io as spio
from optimize import optimizer
from multiprocessing import Pool
import time

        
def procParallel(xT,yT,zT,xB, nProcs=8):
        ''' parallel processing routine using the multiprocess library'''
        pool = Pool(processes=nProcs)
        tic = time.time()
        q = pool.map_async(wrapR5p,zip(xT,yT,zT,xB))
        q.wait()
        print 'Parallel runtime, ' + repr(nProcs) + ' ' + repr(time.time()-tic)
        
        x = np.zeros(xT.size,dtype='complex128')
        y = np.zeros(yT.size,dtype='complex128')
        z = np.zeros(zT.size,dtype='complex128')
        
        for ixr,res in enumerate(q.get()):
            x[ixr] = res[0]
            y[ixr] = res[1]
            z[ixr] = res[2]
            
        pool.terminate()
        return x,y,z
        
        
def wrapR5p(T):
    return r5p(T[0],T[1],T[2],T[3])

    
def r5p(xT,yT,zT,xB):
    '''routine to actually compute the projection '''
    sttm = time.time()
    A0 = np.eye(5)
    Ar = np.zeros((5,5))
    Ai = np.zeros((5,5))
    Ar[0,4] = 0.5
    Ar[4,0] = 0.5
    Ai[1,4] = 0.5
    Ai[4,1] = 0.5
    b0 = -np.array([xT.real, xT.imag, zT.real, zT.imag,yT.real])
    br = np.zeros(5) 
    br[2] = -0.5 
    br[4] = 0.5*xB.real
    bi = np.zeros(5)
    bi[3] = -0.5
    bi[4] = 0.5*xB.imag
    
    F0 = -np.vstack((np.hstack((A0,b0.reshape(5,1))),np.hstack((b0,0.0))))
    F1 = -np.vstack((np.hstack((Ar,br.reshape(5,1))),np.hstack((br,0.0))))
    F2 = -np.vstack((np.hstack((Ai,bi.reshape(5,1))),np.hstack((bi,0.0))))
    F3 = np.zeros((6,6))
    F3[5,5] = -1.0
    
    FM = np.hstack((F1.reshape(36,1), F2.reshape(36,1),F3.reshape(36,1)))
    
    ''' initialize some variables '''
    t = 1.0
    L = np.zeros(3)
    L[2] = 1.0
    c = np.array([0.0, 0.0, 1.0])
    
    ''' find a feasible L -- you can always find one (usually)'''
    stp = 0
    S = -F0 - F1*L[0] - F2*L[1] -F3*L[2]
    while np.any(np.linalg.eig(S)[0] < 0 ) & (stp<10):
        L *= 2
        S = -F0 - F1*L[0] - F2*L[1] -F3*L[2]
        stp += 1
    
#        print 'S ' + repr(S)
    g = np.zeros(3)
    P = []
    H = np.zeros((3,3))
    for itr in range(20):
        SI = np.linalg.pinv(S)
#            print 'SI ' + repr(SI)
        P = ['','','']
        for ixi in range(3):
#                print repr(ixi) + ' at ' + repr(FM[:,ixi].reshape(6,6))
            
            P[ixi] = np.dot(SI,FM[:,ixi].reshape(6,6))
            g[ixi] = c[ixi]*t + np.trace(P[ixi])
            H[ixi,ixi] = np.trace(np.dot(P[ixi],P[ixi]))
#                print H[ixi,ixi]
#                print repr(ixi) + ' P ' + repr(P[ixi])
            
        for ixi in range(3):
            for ixj in range(ixi+1,3):
#                    print repr(ixi) + ' ' + repr(ixj)
                H[ixi,ixj] = np.trace(np.dot(P[ixi],P[ixj]))
                H[ixj,ixi] = np.trace(np.dot(P[ixj],P[ixi]))
#                    print H
        
#            print P
        HI = np.linalg.pinv(H)
        dy = np.dot(HI,-g)

        aa = 1.0
        btk = 0
        go = True
        while go & (btk<15):
            Ln = L + aa*dy
#                print Ln.shape
            Sn = -F0 - F1*Ln[0] - F2*Ln[1] - F3*Ln[2]
            cng = np.linalg.norm(Ln-L)/np.linalg.norm(L)
            if np.any(np.linalg.eig(Sn)[0]<0):
                aa = aa/2
            else:
                L = Ln
                break
            
            btk += 1
#            print btk
        if cng < 1e-6:
            break
        
        t = t*2.5
        S = -F0 - F1*L[0] - F2*L[1] - F3*L[2]
        ''' end of interior point loop'''
    
    M = A0 + L[0]*Ar + L[1]*Ai
    b = b0 + L[0]*br + L[1]*bi
    M = np.linalg.pinv(M)
    
    xh = -np.dot(M,b)
    x = xh[0] + 1j*xh[1]
    y = xh[4]
    z = xh[2] + 1j*xh[3]
    gap = np.linalg.norm((x+xB)*y-z)
    print 'RT = ' + repr(time.time()-sttm) + ' Gap ' + repr(gap)
    return x,y,z
    
        
class problem(optimizer):
    ''' class that extents the contrast - Xadmm algorithm '''
    def initOpt(self, uHat, D):
        self.rho = D['rho']
        self.xi = D['xi']
        self.uHat = uHat
        self.upperBound = D['uBound']
        self.lmb = D['lmb']
        self.obj = np.zeros(D['maxIter'])
        
        # add some local vars for ease
        self.s = self.fwd.getS() #  1j*self.muo*self.w
        self.A = self.fwd.nabla2+self.fwd.getk(0)
        
        # variables for the scattered fields
        self.us = np.zeros(self.fwd.N,dtype='complex128')
        self.uT = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.uD = np.zeros(self.fwd.getXSize(),dtype='complex128')
        
        # variables for the contrast source
        self.x = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.xT = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.xD = np.zeros(self.fwd.getXSize(),dtype='complex128')
        
        # variables for theta
        self.tT = np.zeros(self.fwd.nRx*self.fwd.nRy, dtype='complex128')
        self.tD = np.zeros(self.fwd.nRx*self.fwd.nRy, dtype='complex128')
        
        
        # just to make life easier:
        self.ub = self.fwd.sol[0] # shouldn't need --> .flatten()
        self.uHat = uHat - self.fwd.Ms*self.fwd.sol[0]
        # create some new operators for doing what is necessary for the 
        # contrast X work
        self.fwd.setCTRX()
        
        # self.indefinite = projector()
        uu = self.fwd.Ms.T*self.fwd.Ms + self.fwd.Md.T*self.fwd.Md*self.rho
        ux = sparse.coo_matrix((self.fwd.N,self.fwd.getXSize()))
        ul = self.A.T.conj()
        
        xx = sparse.eye(self.fwd.getXSize(),self.fwd.getXSize())*self.rho
        xl = self.fwd.x2u.T
        
        ll = sparse.coo_matrix((self.fwd.N, self.fwd.N))
                
        M = spt.vCat([spt.hCat([uu, ux, ul]), \
                      spt.hCat([ux.T, xx, xl]),\
                      spt.hCat([ul.T.conj(), xl.T.conj(), ll])])
        
        print M.shape
        self.aux = lin.factorized(M)
        
        
    
    def runOpt(self,P):
        ''' local set of iterations '''
        ''' P updated ahead of here, except in first iteration '''
        ''' update u,x '''
        rhs = np.concatenate((self.fwd.Ms.T*self.uHat + self.rho*self.fwd.x2u*(self.uT - self.uD),\
                         self.rho*(self.xT - self.xD),\
                         np.zeros(self.fwd.N)))
        
        updt = self.aux(rhs)
        
        self.us = updt[:self.fwd.N]
        self.x = updt[self.fwd.N:(self.fwd.N+self.fwd.getXSize())]
        
        ''' update uT,xT,tT '''
        uL = self.fwd.x2u.T*self.us
        uLb = self.fwd.x2u.T*self.ub
#        nn = self.fwd.getXSize()
        
        
        self.uT,self.tT,self.xT = procParallel(uL+self.uD, P+self.tD, (self.x+self.xD)/self.s, uLb)
        self.xT = self.xT*self.s
#        
#        for ix in range(nn):
#            self.uT[ix],self.tT[ix],self.xT[ix] = r5p(uL[ix]+self.uD[ix],\
#                                                                      P[ix]+ self.tD[ix],\
#                                                                      (self.x[ix]+self.xD[ix])/self.s,\
#                                                                      uLb[ix])
#        
        
        ''' update the dual vars '''
        self.uD += uL - self.uT
        self.tD += P - self.tT
        self.xD += self.x - self.xT
        
        obj = np.linalg.norm(self.uHat-self.fwd.Ms*self.us)
        return obj
        
    def writeOut(self, rank, ix=0):
        import os
        assert os.path.exists(self.outDir + 'Data')
        
        us = self.fwd.parseFields(self.us)
        ub = self.fwd.parseFields(self.fwd.sol[0])
        sgmm = self.fwd.parseFields(self.fwd.sigmap[0])
        uTrue = self.fwd.parseFields(self.fwd.sol[1])
            
        D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgmm[0], 'ub':ub[0], \
             'us':us[0], 'uTrue':uTrue[0], 'uT':self.uT, \
             'x':self.x, 'xT':self.xT, 'tD':self.tD,'tT':self.tT, 'obj':self.obj, 'flavor':self.fwd.flavor}
        
        spio.savemat(self.outDir + 'Data/contrastX' + repr(rank) + '_' + repr(ix), D)
        

#    
    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! '''
        N = np.size(S)
        nX = self.fwd.getXSize()
        tL = np.zeros(nX,dtype='complex128')
        
        # we just have to do basic averaging
        for L in S:
            tL += L.tT - L.tD
        
        tL = tL/N
        T = np.zeros(nX,dtype='complex128')
        
        T = comm.allreduce(tL,T,op=MPI.SUM)
        T = T/comm.Get_size() # take the mean!
        T = T.real
        T = np.maximum(T,0)
        # print self.upperBound
        T = np.minimum(T,self.upperBound)
        # print P[1]
               
        return T
        
    def plotSemiParallel(self,P,resid,rank,ix=0):
        ''' Plotting routine if things are semiParallel'''
        import matplotlib.pyplot as plt
        plt.close('all')
        import os
        
        plt.close('all')
        if not os.path.exists(self.outDir + 'Figs'):
            os.mkdir(self.outDir + 'Figs')
        
        # vv = S.Ms*S.v
        uu = self.fwd.Ms*self.us
        ub = self.fwd.Ms*self.ub
        skt = self.uHat-ub
        
        plt.figure(100 + rank + 10*ix)
        plt.plot(np.arange(self.fwd.nSen), skt.real, np.arange(self.fwd.nSen), uu.real)
        plt.savefig(self.outDir + 'Figs/fig' + repr(100+rank+10*ix))
        
        if (rank==0) & (ix==0):
            # then print some figures   
            plt.figure(383)
            plt.plot(resid)
            plt.savefig(self.outDir + 'Figs/fig383')
        
            plt.figure(387)
            plt.imshow(P.reshape(self.fwd.nRx,self.fwd.nRy), interpolation='nearest')
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig387')
            
            us = self.fwd.parseFields(self.us)
            
            plt.figure(76)
            plt.subplot(121)
            plt.imshow(us[0].real)
            plt.colorbar()
        
            plt.subplot(122)
            plt.imshow(us[0].imag)
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig76')
        
        # all show!
        # plt.show()
