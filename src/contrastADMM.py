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
from superSolve import wrapCvxopt
import time

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
        
        self.gap = list()
        self.objInt = list()
        self.pL = list()
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
        self.scaleC = 1.0 # 1.0/np.linalg.norm(self.fwd.Ms*self.ub)
        # create some new operators for doing what is necessary for the 
        # contrast X work
        self.X = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.Z = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.fwd.setCTRX()
        
        ''' subtract out the background field '''
        self.uHat = self.uHat - self.fwd.Ms*self.ub
        # pfake = (self.upperBound/2.0)*np.ones(self.fwd.getPSize(),dtype='complex128')
        ''' in this instance, I don't care about the results, i.e., I don't care about the actual solutions'''
        # self.internalSymbolic(pfake)
    
    def internalSymbolic(self,thk):
        '''create an internal method that 
        (1) knows the structure of the matrix
        (2) only needs new theta estimates
        (3) keeps the symbolic factorization to reuse '''
        cttm = time.time()
        nX = self.fwd.getXSize()
        pm = sparse.spdiags(self.s*self.fwd.p2x*thk, 0, nX, nX)
        # print pm.shape
        # print self.fwd.x2u.shape
        ''' ds changes at ever iteration '''
        ds = pm*self.fwd.x2u.T #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        ''' changes '''
        bmuu = self.scaleC*(self.fwd.Ms.T*self.fwd.Ms) + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmxu = -self.rho*ds 
        ''' static ''' 
        bmul = self.A.T.conj()
        bmxx = self.rho*sparse.eye(nX, nX)
        bmxl = self.fwd.x2u.T
        bmlu = self.A
        bmlx = self.fwd.x2u
        bmll = sparse.coo_matrix((self.fwd.N, self.fwd.N))
        
        ''' right hand side ''' 
        rhsu = self.scaleC*self.fwd.Ms.T.conj()*self.uHat - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        rhsx = self.rho*ds*self.ub - self.rho*self.Z # chng
        rhsl = np.zeros(self.fwd.N)
  
  
        bm = spTools.vCat([spTools.hCat([bmuu, bmux, bmul]), \
                           spTools.hCat([bmxu, bmxx, bmxl]), \
                           spTools.hCat([bmlu, bmlx, bmll])])
        
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        print 'construction time ' + repr(time.time()-cttm)
        
        sltm = time.time()
        if hasattr(self,'symbFact'):
            print 'Ive already got it'
            updt = wrapCvxopt.solveNumeric(bm, rhsbm, self.symbFact)
        else:
            print 'got to do the symbolic calc'
            self.symbFact = wrapCvxopt.createSymbolic(bm)
            updt = wrapCvxopt.solveNumeric(bm, rhsbm, self.symbFact)
        
        # updt = lin.spsolve(bm.tocsr(), rhsbm)
        
        # N = self.nx*self.ny
        us = updt[:self.fwd.N]
        x = updt[self.fwd.N:(self.fwd.N+nX)]
        print 'solve time ' + repr(time.time()-sltm)
        return us,x
        
    def internalHard(self, thk):
        '''creates the matrix every time, a hard alternative, to internalSymbolic
        so that I can do A/B testing easily w.r.t the old standard'''
        nX = self.fwd.getXSize()
        pm = sparse.spdiags(self.s*self.fwd.p2x*thk, 0, nX, nX)
        # print pm.shape
        # print self.fwd.x2u.shape
        ''' ds changes at ever iteration '''
        ds = pm*self.fwd.x2u.T #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        ''' changes '''
        bmuu = (self.fwd.Ms.T*self.fwd.Ms) + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmxu = -self.rho*ds 
        ''' static ''' 
        bmul = self.A.T.conj()
        bmxx = self.rho*sparse.eye(nX, nX)
        bmxl = self.fwd.x2u.T
        bmlu = self.A
        bmlx = self.fwd.x2u
        bmll = sparse.coo_matrix((self.fwd.N, self.fwd.N))
        
        ''' right hand side ''' 
        rhsu = self.fwd.Ms.T.conj()*self.uHat - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        rhsx = self.rho*ds*self.ub - self.rho*self.Z # chng
        rhsl = np.zeros(self.fwd.N)
  
  
        bm = spTools.vCat([spTools.hCat([bmuu, bmux, bmul]), \
                           spTools.hCat([bmxu, bmxx, bmxl]), \
                           spTools.hCat([bmlu, bmlx, bmll])])
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        
        
        
        updt = lin.spsolve(bm.tocsr(), rhsbm)
        if self.fwd.flavor == 'TE3D':
            spio.savemat('hugeM', {'bm':bm, 'rhs':rhsbm, 'sol':updt})
        elif self.fwd.flavor == 'TE':
            spio.savemat('normalTE', {'bm':bm, 'rhs':rhsbm, 'sol':updt})
            
        
        # N = self.nx*self.ny
        us = updt[:self.fwd.N]
        x = updt[self.fwd.N:(self.fwd.N+nX)]
        return us,x
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        
        '''update dual variables first '''
        print 'size of us ' + repr(self.us.shape)
        print 'size of x2u ' + repr(self.fwd.x2u.shape)
        print 'size of p2x ' + repr(self.fwd.p2x.shape)
        print 'size of Z ' + repr(self.Z.shape) + ' and X ' + repr(self.X.shape)
        self.Z = self.Z + (self.X - (self.s*self.fwd.x2u.T*(self.ub + self.us))*(self.fwd.p2x*P))
        
        ''' jointly update u,x '''
        # pfake = (self.upperBound/2.0)*np.ones(self.fwd.getXSize(),dtype='complex128')
        self.us,self.X = self.internalHard(P)
        
        self.gap.append(np.linalg.norm(self.X - (self.fwd.p2x*P)*(self.s*self.fwd.x2u.T*(self.us+self.ub))))
        
        obj = np.linalg.norm(self.uHat-self.fwd.Ms*self.us)
        self.objInt.append(obj)
        return obj
        
    def writeOut(self, rank, ix=0):
        import os
        assert os.path.exists(self.outDir + 'Data')
        
        us = self.fwd.parseFields(self.us)
        ub = self.fwd.parseFields(self.fwd.sol[0])
        sgmm = self.fwd.parseFields(self.fwd.sigmap[1])
        uTrue = self.fwd.parseFields(self.fwd.sol[1])
        
        if (self.fwd.flavor == 'TE') | (self.fwd.flavor == 'TM'):   
            D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgmm[0], 'ub':ub[0], \
                 'us':us[0], 'uTrue':uTrue[0], \
                 'X':self.X, 'obj':self.obj, 'flavor':self.fwd.flavor, 'gap':self.gap, \
                 'obj':self.objInt, 'Ms':self.fwd.Ms, 'phist':self.pL}
        elif (self.fwd.flavor == 'TE3D'):
            D = {'f':self.fwd.f, 'angle':self.fwd.incAng, \
                 'sigMatX':sgmm[0], 'sigMatY':sgmm[1], 'sigMatZ': sgmm[2],\
                 'ubX':ub[0], 'ubY':ub[1], 'ubZ':ub[2], \
                 'usX':us[0], 'usY':us[1], 'usZ':us[2], \
                 'uTrueX':uTrue[0], 'uTrueY':uTrue[1], 'uTrueZ':uTrue[2], \
                 'X':self.X, 'obj':self.obj, 'flavor':self.fwd.flavor, 'gap':self.gap, \
                 'obj':self.objInt, 'Ms':self.fwd.Ms, 'phist':self.pL,'Md':self.fwd.Md}
        
        spio.savemat(self.outDir + 'Data/contrastX' + repr(rank) + '_' + repr(ix), D)
        
        

    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! '''
        ''' oh this is going to get trickier when we return to 3d & multivariates'''
        n = S[0].fwd.getPSize()
        nX = S[0].fwd.getXSize()
        
        ''' matrix method, i.e. TM, TE3D req'd' '''
        uL = sparse.lil_matrix((n,n),dtype='complex128')
        bL = np.zeros(n,dtype='complex128')
        
        for L in S:
            M = L.s*(sparse.spdiags(L.fwd.x2u.T*(L.ub+L.us),0,nX,nX))*self.fwd.p2x
            uL = uL + M.T.conj()*M
            bL = bL + M.T.conj()*(L.X + L.Z)
            
        U = sparse.lil_matrix((n,n),dtype='complex128')
        B = np.zeros(n,dtype='complex128')

        U = comm.allreduce(uL,U,op=MPI.SUM)
        B = comm.allreduce(bL,B,op=MPI.SUM)
        ''' interesting: += isn't implemented?'''
        U = U + sparse.spdiags((self.lmb/self.rho)*np.ones(n), 0, n, n)
        
        P = lin.spsolve(U,B)
        self.pL.append(lin.spsolve(uL,bL))
        
        ''' matrix free method '''
#        uL = np.zeros(n,dtype='complex128')
#        bL = np.zeros(n,dtype='complex128')
#        for L in S:
#            uL += self.s*self.fwd.Md*(self.us + self.ub)
#            bL += self.X + self.Z

#        U = np.zeros(n,dtype='complex128')
#        B = np.zeros(n,dtype='complex128')
#        
#        U = comm.allreduce(uL,U,op=MPI.SUM)
#        B = comm.allreduce(bL,B,op=MPI.SUM)
#        
#        num = self.rho*(U.real*B.real + U.imag*B.imag)
#        den = self.rho*(U.conj()*U) + self.lmb
#        
#        P = (num/den)
#        self.pL.append((uL.real*bL.real + uL.imag*bL.imag)/(uL.conj()*uL))        
        
        ''' finish off '''
        P = P.real
        P = np.maximum(P,0)
        P = np.minimum(P,self.upperBound)
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
