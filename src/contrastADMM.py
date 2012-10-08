'''
Created on Jun 3, 2012

@author: dstrauss

implementation of contrast source ADMM optimization
'''

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import sparseTools as spTools
import scipy.io as spio
from optimize import optimizer

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
        
        # oh. this is going to get strange. 
        # integrate TE concepts first.
        # contrast variable
#        self.X = np.zeros(self.fwd.nRx*self.fwd.nRy,dtype='complex128')
        # dual variable 
#        self.Z = np.zeros(self.fwd.nRx*self.fwd.nRy,dtype='complex128')
        # scattered fields
        self.us = np.zeros(self.fwd.N,dtype='complex128')
        # just to make life easier:
        self.ub = self.fwd.sol[0] # shouldn't need --> .flatten()
        
        # create some new operators for doing what is necessary for the 
        # contrast X work
        self.X = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.Z = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.fwd.setCTRX()
        
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        self.Z = self.Z + (self.X - (self.s*self.fwd.x2u.T*(self.ub + self.us))*(self.fwd.p2x*P))
        
        uHatLocal =  self.uHat - self.fwd.Ms*self.ub  #remove background field
        
        nX = self.fwd.getXSize()
        pm = sparse.spdiags(self.s*self.fwd.p2x*P, 0, nX, nX)
        # print pm.shape
        # print self.fwd.x2u.shape
        ds = pm*self.fwd.x2u.T #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        bmuu = self.fwd.Ms.T*self.fwd.Ms + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmul = self.A.T.conj()
        
        rhsu = self.fwd.Ms.T.conj()*uHatLocal - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        
  
        bmxu = -self.rho*ds
        bmxx = self.rho*sparse.eye(nX, nX)
        bmxl = self.fwd.x2u.T
        rhsx = self.rho*ds*self.ub - self.rho*self.Z
  
        bmlu = self.A
        bmlx = self.fwd.x2u

        bmll = sparse.coo_matrix((self.fwd.N, self.fwd.N)) 
        rhsl = np.zeros(self.fwd.N)
  
  
        bm = spTools.vCat([spTools.hCat([bmuu, bmux, bmul]), \
                             spTools.hCat([bmxu, bmxx, bmxl]), \
                             spTools.hCat([bmlu, bmlx, bmll])])
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        
        updt = lin.spsolve(bm.tocsr(), rhsbm)
        
        # N = self.nx*self.ny
        self.us = updt[:self.fwd.N]
        self.X = updt[self.fwd.N:(self.fwd.N+nX)]
        
        obj = np.linalg.norm(uHatLocal-self.fwd.Ms*self.us)
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
        n = S[0].fwd.nRx*S[0].fwd.nRy
        
        uL = sparse.lil_matrix((n,n),dtype='complex128')
        bL = np.zeros(n,dtype='complex128')
        
#        U = np.zeros((n,N),dtype='complex128')
#        Q = np.zeros((n,N),dtype='complex128')
        nX = self.fwd.getXSize()
        for L in S:
            # for ix in range(N):
#            s = S[ix].s
            # print s
#            U[:,ix] = s*S[ix].fwd.Md*(S[ix].ub + S[ix].us)
#            Q[:,ix] = S[ix].X + S[ix].Z
            # print L.fwd.x2u.shape
            
            M = L.s*(sparse.spdiags(L.fwd.x2u.T*(L.ub+L.us),0,nX,nX))*self.fwd.p2x
            uL += M.T.conj()*M
            bL += M.T.conj()*(L.X + L.Z)
            
            
#        numLocal = np.sum(U.real*Q.real + U.imag*Q.imag,1)
#        denLocal = np.sum(U.conj()*U,1) + self.lmb/S[0].rho
#        
#        num = np.zeros(numLocal.shape)
#        num = comm.allreduce(numLocal,num,op=MPI.SUM)
#
#        den = np.zeros(denLocal.shape)
#        den = comm.allreduce(denLocal,den,op=MPI.SUM)
        U = sparse.lil_matrix((n,n),dtype='complex128')
        B = np.zeros(n,dtype='complex128')
        
        U = comm.allreduce(uL,U,op=MPI.SUM)
        B = comm.allreduce(bL,B,op=MPI.SUM)
        
        P = lin.spsolve(U,B)
        # print num
        # print den
        
        P = P.real
        
        # print P[1]
        P = np.maximum(P,0)
        # print self.upperBound
        P = np.minimum(P,self.upperBound)
        # print P[1]
               
        return P
        
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
