'''
Created on Jun 6, 2012

@author: dstrauss
'''

import numpy as np
from optimize import optimizer
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
import scipy.io as spio
from mpi4py import MPI

class problem(optimizer):
    '''Implementing vanilla biconvex method for subsurface imaging '''
    def initOpt(self, uHat, rho, xi, upperBound, lmb, maxiter):
        '''initialization routine with common prototype'''
        self.rho = rho
        self.xi = xi
        self.uHat = uHat
        self.lmb = lmb
        self.uBound = upperBound
        self.F = np.zeros(self.fwd.N, dtype='complex128')
        self.us = np.zeros(self.fwd.N, dtype='complex128')
        
        self.A = self.fwd.nabla2 + self.fwd.getk(0)
        self.s = self.fwd.getS()
        self.ub = self.fwd.sol[0]
        self.obj = np.zeros(maxiter)
        
        # not much else to do but initialize variables
        
    def runOpt(self,P):
        ''' inner loop protocol '''
        
        # update dual variables
        self.F += (self.A*self.us + (self.s*self.fwd.Md.T*P)*(self.ub + self.us))
        uHatLocal = self.uHat - self.fwd.Ms*self.ub
        aL  = self.A + sparse.spdiags(self.s*self.fwd.Md.T*P, 0, self.fwd.N, self.fwd.N)
        M = self.fwd.Ms.T*self.fwd.Ms + self.rho*aL.T.conj()*aL
        
        b = self.fwd.Ms.T*uHatLocal - \
            self.rho*(aL.T.conj()*(self.s*self.ub*(self.fwd.Md.T*P) + self.F ))
        
        # solve dat dere set of equations
        self.us = lin.spsolve(M,b)
        
        obj = np.linalg.norm(self.fwd.Ms*self.us - uHatLocal)
        return obj
    
    def writeOut(self,rank,ix=0):
        
        ub = self.fwd.parseFields(self.fwd.sol[0])
        ut = self.fwd.parseFields(self.fwd.sol[1])
        
        us = self.fwd.parseFields(self.us)
        
        D = {'f':self.fwd.f, 'angle': self.fwd.incAng, 'ub':ub[0], \
             'us':us[0], 'uTrue':ut[0], 'rho':self.rho, 'xi':self.xi, 'obj':self.obj}
        spio.savemat(self.outDir + 'Data/biconvex' + repr(rank) + '_' + repr(ix), D)
        
    def aggregatorSemiParallel(self,S,comm):
        '''semiParallel method does aggregation over the local list, then distributes
        through the cluster'''
        
        N = np.size(S)
        n = self.fwd.nRx*self.fwd.nRy
        
        uL = np.zeros(n,dtype='complex128')
        qL = np.zeros(n,dtype='complex128')
        
        for L in S:
            uL += L.s*L.fwd.Md*(self.ub + self.us)
            qL += L.fwd.Md*(L.A*(self.us + self.F))

        uL = uL*(1.0/N)
        qL = qL*(1.0/N)
        
        U = np.zeros(n,dtype='complex128')
        Q = np.zeros(n,dtype='complex128')
        
        U = comm.allreduce(uL,U,op=MPI.SUM)
        Q = comm.allreduce(qL,Q,op=MPI.SUM)
        
        U = U*(1.0/comm.Get_size())
        Q = Q*(1.0/comm.Get_size())
        
        num = self.rho*(U.real*Q.real + U.imag*Q.imag)
        den = self.rho*(U.conj()*U) + self.lmb

        P = (-num/den).real
        P = np.maximum(P,0)
        P = np.minimum(P,self.uBound)
        
        return P
    
    def plotSemiParallel(self, P, resid, rank,ix=0):
        ''' plotting routine for the semiParallel case -- this should probably just be made
        standard, esp. since it is rarely used now '''
        import matplotlib.pyplot as plt
        plt.close('all')
        
        uu = self.fwd.Ms*self.us
        ub = self.fwd.Ms*self.ub
        skt = self.uHat-ub
        
        plt.figure(100+rank+10*ix)
        plt.plot(np.arange(self.fwd.nSen), skt.real, np.arange(self.fwd.nSen), uu.real)
        plt.savefig(self.outDir + 'Figs/fig' + repr(100+rank+10*ix) )
        
        if rank==0 & ix==0:
            # then print some figures   
            plt.figure(383)
            plt.plot(resid)
            plt.savefig(self.outDir + 'Figs/fig383'  )
        
            plt.figure(387)
            plt.imshow(P.reshape(self.fwd.nRx, self.fwd.nRy), interpolation='nearest')
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig387'  )
            
            u = self.fwd.parseFields(self.us)
            plt.figure(76)
            plt.subplot(121)
            plt.imshow(u[0].real)
            plt.colorbar()

            plt.subplot(122)
            plt.imshow(u[0].imag)
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig76'  )
        
        
        