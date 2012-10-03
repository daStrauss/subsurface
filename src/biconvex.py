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
# import superSolve.wrapCvxopt

class problem(optimizer):
    '''Implementing vanilla biconvex method for subsurface imaging '''
    def initOpt(self, uHat, D):
        '''initialization routine with common prototype'''
        self.rho = D['rho']
        self.xi = D['xi']
        self.uHat = uHat
        self.lmb = D['lmb']
        self.uBound = D['uBound']
        self.F = np.zeros(self.fwd.N, dtype='complex128')
        self.us = np.zeros(self.fwd.N, dtype='complex128')
        
        self.A = self.fwd.nabla2 + self.fwd.getk(0)
        self.s = self.fwd.getS()
        self.ub = self.fwd.sol[0]
        self.obj = np.zeros(D['maxIter'])
        
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
#        self.us = superSolve.wrapCvxopt.linsolve(M, b)
        
        obj = np.linalg.norm(self.fwd.Ms*self.us - uHatLocal)
        return obj
    
    def writeOut(self,rank,ix=0):
        
        ub = self.fwd.parseFields(self.fwd.sol[0])
        ut = self.fwd.parseFields(self.fwd.sol[1])
        
        us = self.fwd.parseFields(self.us)
        
        D = {'f':self.fwd.f, 'angle': self.fwd.incAng, 'ub':ub[0], \
             'us':us[0], 'uTrue':ut[0], 'rho':self.rho, 'xi':self.xi, 'obj':self.obj,\
             'flavor':self.flavor}
        spio.savemat(self.outDir + 'Data/biconvex' + repr(rank) + '_' + repr(ix), D)
        
    def aggregatorSemiParallel(self,S,comm):
        '''semiParallel method does aggregation over the local list, then distributes
        through the cluster'''
        
        N = np.size(S)
        n = self.fwd.nRx*self.fwd.nRy # ie size of P.
        
        uL = np.zeros(n,dtype='complex128')
        qL = np.zeros(n,dtype='complex128')

#        aL = sparse.lil_matrix((n,n),dtype='complex128')
#        bL = np.zeros(n,dtype='complex128')
        
        for L in S:
#            D = sparse.spdiags(L.s*(L.ub+L.us),0, self.fwd.N,self.fwd.N)
#            g = D*L.fwd.Md.T
#            bL += g.T.conj()*((L.A*L.us) + L.F)
#            aL += g.T.conj()*g
            
             uL += L.s*L.fwd.Md*(self.ub + self.us)
             qL += L.fwd.Md*(L.A*self.us + self.F)

        # uL = uL*(1.0/N)
        # qL = qL*(1.0/N)
        
#        A = sparse.lil_matrix((n,n),dtype='complex128')
#        B = np.zeros(n,dtype='complex128')
        
        U = np.zeros(n,dtype='complex128')
        Q = np.zeros(n,dtype='complex128')
        
        U = comm.allreduce(uL,U,op=MPI.SUM)
        Q = comm.allreduce(qL,Q,op=MPI.SUM)
#        A = comm.allreduce(aL,A,op=MPI.SUM)
#        B = comm.allreduce(bL,B,op=MPI.SUM)
        # U = U*(1.0/comm.Get_size())
        # Q = Q*(1.0/comm.Get_size())
        
        num = self.rho*(U.real*Q.real + U.imag*Q.imag)
        den = self.rho*(U.conj()*U) + self.lmb
#        A += sparse.eye(n,n)*self.lmb/self.rho
        
#        P = lin.spsolve(A.tocsr(),-B).real
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
        
        
        