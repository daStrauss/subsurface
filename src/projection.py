'''
Created on Oct 22, 2012

@author: dstrauss

Implementation of the projection operator method for solving the subsurface imaging problem
'''

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import sparseTools as spt
import scipy.io as spio
from optimize import optimizer

class projector(object):
    ''' class for implementing the R5 projection over (x+xb)y=z ''' 
   
    
    def __init__(self):
        '''ok make some internal vars that will get reused as time goes along '''
        self.A0 = np.eye(5)
        self.Ar = np.zeros((5,5))
        self.Ai = np.zeros((5,5))
        self.Ar[0,4] = 0.5
        self.Ar[4,0] = 0.5
        self.Ai[1,4] = 0.5
        self.Ai[4,1] = 0.5
    
    def r5p(self,xT,yT,zT,xB):
        '''routine to actually compute the projection '''
        b0 = -np.array([xT.real, xT.imag, zT.real, zT.imag,yT.real])
        br = np.zeros(5) 
        br[2] = -0.5 
        br[4] = 0.5*xB.real
        bi = np.zeros(5)
        bi[3] = -0.5
        bi[4] = 0.5*xB.imag
        
        F0 = -np.vstack((np.hstack((self.A0,b0.reshape(5,1))),np.hstack((b0,0.0))))
        F1 = -np.vstack((np.hstack((self.Ar,br.reshape(5,1))),np.hstack((br,0.0))))
        F2 = -np.vstack((np.hstack((self.Ai,bi.reshape(5,1))),np.hstack((bi,0.0))))
        F3 = np.zeros((6,6))
        F3[5,5] = -1.0
        
#        print F0
#        print F1
#        print F2
#        print F3
        
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
        
        M = self.A0 + L[0]*self.Ar + L[1]*self.Ai
        b = b0 + L[0]*br + L[1]*bi
        M = np.linalg.pinv(M)
        
        xh = -np.dot(M,b)
        x = xh[0] + 1j*xh[1]
        y = xh[4]
        z = xh[2] + 1j*xh[3]
        gap = np.linalg.norm((x+xB)*y-z)
        print 'Gap ' + repr(gap)
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
        
        self.indefinite = projector()
        uu = self.fwd.Ms.T*self.fwd.Ms + self.fwd.Md.T*self.fwd.Md*self.rho
        ux = sparse.coo_matrix((self.fwd.N,self.fwd.getXSize()))
        ul = self.A.T.conj()
        
        xx = sparse.eye(self.fwd.getXSize(),self.fwd.getXSize())*self.rho
        xl = self.fwd.x2u.T
        
        ll = sparse.coo_matrix(self.fwd.N, self.fwd.N)
        
        print uu.shape
        print ux.shape
        print ul.shape
        
        print xx.shape
        print xl.shape
        
        print ll.shape
        
        M = spt.vCat([spt.hCat([uu, ux, ul]), \
                      spt.hCat([ux.T, xx, xl]),\
                      spt.hCat([ul.T.conj(), xl.T.conj, ll])])
        
        self.aux = lin.factorized(M)
        
        
    
    def runOpt(self,P):
        ''' local set of iterations '''
        ''' P updated ahead of here, except in first iteration '''
        ''' update u,x '''
        rhs = np.vstack((self.fwd.Ms.T*self.uHat + self.rho*self.fwd.x2u*(self.uT - self.uD),\
                         self.rho*(self.xT - self.xD),\
                         np.zeros(self.fwd.N)))
        
        updt = self.aux(rhs)
        
        self.us = updt[:self.fwd.N]
        self.x = updt[self.fwd.N:(self.fwd.N+self.fwd.getXSize())]
        
        ''' update uT,xT,tT '''
        uL = self.fwd.x2u.T*self.us
        uLb = self.fwd.x2u.T*self.ub
        nn = self.fwd.getXSize()
        
        for ix in range(nn):
            self.uT[ix],self.tT[ix],self.xT[ix] = self.indefinite.r5p(uL[ix]+self.uD[ix],\
                                                                      P[ix]+ self.tD[ix],\
                                                                      self.x[ix]+self.xD[ix],\
                                                                      uLb[ix])
        
        
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
             'us':us[0], 'uTrue':uTrue[0], \
             'X':self.X, 'obj':self.obj, 'flavor':self.fwd.flavor}
        
        spio.savemat(self.outDir + 'Data/contrastX' + repr(rank) + '_' + repr(ix), D)
        
        
#   Just break it -- I don't plan to use it anyway     
#    def aggregatorSerial(self, S):
#        ''' routine to do the aggregation step and return an update for P '''
#        N = np.size(S)
#        n = S[0].nRx*S[0].nRy
#        # print N
#        # print n
#        
#        U = np.zeros((n,N),dtype='complex128')
#        Q = np.zeros((n,N),dtype='complex128')
#        
#        for ix in range(N):
#            s = S[ix].s
#            U[:,ix] = s*S[ix].Md*(S[ix].ub + S[ix].us)
#            Q[:,ix] = S[ix].X + S[ix].Z
#            
#        num = np.sum(U.real*Q.real + U.imag*Q.imag,1)
#        
#        den = np.sum(U.conj()*U,1) + self.lmb/S[0].rho
#        
#        P = (num/den).real
#        P = np.maximum(P,0)
#        P = np.minimum(P,self.upperBound)
#        
#        gap = np.zeros(N)
#        for ix in range(N):
#            gap[ix] = np.linalg.norm(U[:,ix]*P - S[ix].X)
#            
#        return P
#    
#    def aggregatorParallel(self, comm):
#        ''' Do the aggregation step in parallel whoop! '''
#        N = np.size(self)
#        print repr(N) + ' == better be 1!'
#        
#        U = self.s*self.Md*(self.ub+self.us)
#        Q = self.X + self.Z
#        
#        q = U.real*Q.real + U.imag*Q.imag
#        num = np.zeros(q.shape)
#        
#        num = comm.allreduce(q,num,op=MPI.SUM)
#        
#        q = U.conj()*U + self.lmb/self.rho
#        den = np.zeros(q.shape)
#        den = comm.allreduce(q,den,op=MPI.SUM)
#        
#        P = (num/den).real
#        P = np.maximum(P,0)
#        P = np.minimum(P,self.upperBound)
#        
#        # gap = np.linalg.norm(U*P - self.X)
#        # print 'Proc ' + repr(comm.Get_rank()) + ' gap = ' + repr(gap)
#        
#        return P
#    
    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! '''
        N = np.size(S)
        nX = self.fwd.getXSize()
        tL = np.zeros(nX)
        
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
        T = np.minimum(T.self.upperBound)
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
