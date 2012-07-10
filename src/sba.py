'''
Created on Jun 4, 2012

@author: dstrauss
'''
# maxwell decomissioned
# from maxwell import twoDim
import scipy.sparse as sparse
# import scipy.sparse.linalg as lin
import sparseTools as spt
import numpy as np
from mpi4py import MPI
import scipy.io as spio
from optimize import optimizer
import superSolve.wrapCvxopt 

class problem(optimizer):
    '''a class to do the born approximation iterations '''
    
    #internal parameters of the solver
    eabs = 1e-4
    erel = 1e-3
    rho = 0.0001
    maxit = 2000
    
    def initOpt(self, uHat, rho=0.005, xi=0.9, uBound=0.05, lmb=0, maxiter=100):
        self.stepSize = rho
        self.alpha = xi
        self.uHat = uHat
        self.s = self.fwd.getS()
        self.uBound = uBound
        self.lmb = lmb
        self.obj = np.zeros(maxiter)
        
        
    def trueSolve(self,P):
        ''' an "internal" method to use for obtaining the current value of us '''
        
        self.A = self.fwd.nabla2 + self.fwd.getk(0) + \
          sparse.spdiags(self.s*self.fwd.Md.T*P,0,self.fwd.N, self.fwd.N)
        
#        self.us = lin.spsolve(self.A,self.fwd.rhs)

        self.us = superSolve.wrapCvxopt.linsolve(self.A, self.fwd.rhs)
    
    def runOpt(self,P):
        ''' runtime module'''
        # get the local background -- requires a new solve!
        self.trueSolve(P)
        localuHat = self.uHat - self.fwd.Ms*self.us
        
        # produce some local matrices
        B = self.s*sparse.spdiags(self.us,0,self.fwd.N,self.fwd.N)*self.fwd.Md.T
        c = np.zeros(self.fwd.N)
        Z = self.fwd.Ms
        
        # grab some internal dimensions -- non ROM for now.
        n = self.fwd.N
        m = self.fwd.nRx*self.fwd.nRy
        
        # initialize some empty variables
        v = np.zeros(n, dtype='complex128') # update to fields (du)
        p = np.zeros(m, dtype='complex128') # estimate 1 for materials
        q = np.zeros(m, dtype='complex128') # estimate 2 for materials
        r = np.zeros(m, dtype='complex128') # dual variable for materials
        
        M = spt.vCat([spt.hCat([Z.T*Z, sparse.coo_matrix((n,m)), self.A.T.conj()]), \
                      spt.hCat([sparse.coo_matrix((m,n)), self.rho*sparse.eye(m,m), B.T.conj()]),\
                      spt.hCat([self.A, B, sparse.coo_matrix((n,n))])]).tocsc()
        
        # print 'M format ' + repr(M.format)
#        Q = lin.factorized(M)
        Q = superSolve.wrapCvxopt.staticSolver(M)
        # Foo = gc.get_objects()
        # print Foo
        # Q = M 
        lowb = -P.copy()
        lowb[P-self.stepSize > 0 ] = -self.stepSize
        
        uppb = self.uBound-P
        uppb[P+self.stepSize < self.uBound] = self.stepSize
        
        okgo = True
        iterk = 0
        
        while okgo:
            iterk += 1
            rhs = np.concatenate((Z.T*localuHat, self.rho*(q-r), c))
            updt = Q(rhs)
            # updt = Q*rhs
            
            v = updt[:n]
            p = updt[n:(n+m)]
            
            qold = q.copy()
            q = (p+r).real
            q[q<=lowb] = lowb[q<=lowb]
            q[q>=uppb] = uppb[q>=uppb]
            
            res = p-q
            ser = self.rho*(q-qold)
            
            epri = np.sqrt(m)*self.eabs + self.erel*max(np.linalg.norm(p), np.linalg.norm(q))
            edual = np.sqrt(m)*self.eabs + self.erel*np.linalg.norm(r)
            
            if (np.linalg.norm(res) <= epri) & (np.linalg.norm(ser) <= edual):
                # print 'At ADMM internal limit at iter ' + repr(iterk) 
                okgo = False
            elif (iterk >= self.maxit):
                okgo = False
                # print 'Hit MAXX iterations'
            else:
                okgo = True
                
            r = r + (p-q)
        
        
        self.deltaP = q
        
        obj = np.linalg.norm(localuHat - self.fwd.Ms*v)
        return obj
        
    
#     decommisioned!
#    def aggregatorSerial(self,S):
#        ''' super simple aggregator for sba method'''
#        N = np.size(S)
#        n = S[0].nRx*S[0].nRy
#        
#        P = np.zeros(n)
#        
#        for ix in range(N):
#            P += S[ix].deltaP
#        
#        return (1.0/N)*P*self.alpha
#    
#    def aggregatorParallel(self,comm):
#        ''' SUPER simple aggregator for sba method for working over parallel '''
#        dP = np.zeros(self.nRx*self.nRy)
#        
#        dP = comm.allreduce(self.deltaP,dP,op=MPI.SUM)
#        dP = dP*(1.0/comm.Get_size())*self.alpha
#        return dP
    
    def aggregatorSemiParallel(self, S, comm):
        ''' super simple aggregator for sba method'''
        # I think that this one stays easy easy
        N = np.size(S)
        n = self.fwd.nRx*self.fwd.nRy
        
        P = np.zeros(n)
        for ix in range(N):
            P += S[ix].deltaP
        
        P = (1.0/N)*P
        dP = np.zeros(self.fwd.nRx*self.fwd.nRy)
        dP = comm.allreduce(P,dP,op=MPI.SUM)
        dP = dP*(1.0/comm.Get_size())*self.alpha
        
        return dP
    
#    def plotParallel(self,P,resid,rank):
#        ''' Plotting routine if things are parallel'''
#        import matplotlib.pyplot as plt
#        import os
#        
#        if not os.path.exists(self.outDir + 'Figs'):
#            os.makedirs(self.outDir + 'Figs')
#            
#        
#        # vv = S.Ms*S.v
#        self.trueSolve(P)
#            
#        uu = self.Ms*(self.us - self.sol[0].flatten())
#        ub = self.Ms*(self.sol[0].flatten())
#        skt = self.uHat-ub
#        
#        plt.figure(100+rank)
#        plt.plot(np.arange(self.nSen), skt.real, np.arange(self.nSen), uu.real)
#        plt.savefig(self.outDir + 'Figs/fig' + repr(100+rank))
#        
#        if rank==0:
#            # then print some figures   
#            plt.figure(383)
#            plt.plot(resid)
#            plt.savefig(self.outDir + 'Figs/fig383')
#        
#            plt.figure(387)
#            plt.imshow(P.reshape(self.nRx,self.nRy), interpolation='nearest')
#            plt.colorbar()
#            plt.savefig(self.outDir + 'Figs/fig387')
#    
#            plt.figure(76)
#            plt.subplot(121)
#            plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).real)
#            plt.colorbar()
#        
#            plt.subplot(122)
#            plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).imag)
#            plt.colorbar()
#            plt.title('Final Scattered Fields f = ' + repr(self.f))
#            plt.savefig(self.outDir + 'Figs/fig76')
#        
#    
#    def plotSerial(self,S,P,resid):
#        ''' plotting routine for the serial passes '''
#        import matplotlib.pyplot as plt
#        import os
#        
#        if not os.path.exists(self.outDir + 'Figs'):
#            os.makedirs(self.outDir + 'Figs')
#        
#        N = np.size(S)
#        plt.figure(383)
#        plt.plot(resid)
#        plt.savefig(self.outDir + 'Figs/fig383')
#        
#        plt.figure(387)
#        plt.imshow(P.reshape(S[0].nRx,S[0].nRy), interpolation='nearest')
#        plt.colorbar()
#        plt.savefig(self.outDir + 'Figs/fig387')
#        
#        for ix in range(N):
#            plt.figure(50+ix)
#            # vv = S[ix].Ms*S[0].v
#            S[ix].trueSolve(P)
#            
#            uu = S[ix].Ms*(self.us - S[ix].sol[0].flatten())
#            ub = S[ix].Ms*(S[0].sol[0].flatten())
#            skt = S[ix].uHat-ub
#        
#            # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})
#    
#            # plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real, np.arange(S[0].nSen), vv.real)
#            plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real)
#            plt.savefig(self.outDir + 'Figs/fig'+repr(50+ix))
#            
#        plt.figure(76)
#        plt.subplot(121)
#        plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).real)
#        plt.colorbar()
#    
#        plt.subplot(122)
#        plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).imag)
#        plt.colorbar()
#        plt.savefig(self.outDir + 'Figs/fig76')
#        plt.title('Final Scattered Fields f = ' + repr(self.f))
#
#        plt.savefig(self.outDir + 'Figs/fig76')
        # plt.show()
    def writeOut(self, rank, ix=0):
        '''routine to print out information about the solve '''
        import os
        if not os.path.exists(self.outDir + 'Data'):
            os.mkdir(self.outDir + 'Data')
            
        
        sgm = self.fwd.parseFields(self.fwd.sigmap[0])
        ub = self.fwd.parseFields(self.fwd.sol[0])
        us = self.fwd.parseFields(self.us)
        uTrue = self.fwd.parseFields(self.fwd.sol[1])
        
        
        D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgm[0], 'ub':ub[0], \
             'us':us[0], 'uTrue':uTrue[0], 'obj':self.obj}
        
        spio.savemat(self.outDir + 'Data/sba' + repr(rank) + '_' + repr(ix), D)
    
    def plotSemiParallel(self,P,resid,rank,ix=0):
        ''' Plotting routine if things are semiParallel'''
        import matplotlib.pyplot as plt
        plt.close('all')
        import os
        
        assert os.path.exists(self.outDir + 'Figs')
        
        # vv = S.Ms*S.v
        self.trueSolve(P)
            
        uu = self.fwd.Ms*(self.us - self.fwd.sol[0])
        ub = self.fwd.Ms*(self.fwd.sol[0])
        skt = self.uHat-ub
        
        plt.figure(100 + rank + ix*10)
        plt.plot(np.arange(self.fwd.nSen), skt.real, np.arange(self.fwd.nSen), uu.real)
        plt.savefig(self.outDir + 'Figs/fig' + repr(100+rank+ix*10))
        
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
            ub = self.fwd.parseFields(self.fwd.sol[0])
            plt.figure(76)
            plt.subplot(121)
            plt.imshow((us[0]-ub[0]).real)
            plt.colorbar()
        
            plt.subplot(122)
            plt.imshow((us[0] - ub[0]).imag)
            plt.colorbar()
            plt.title('Final Scattered Fields f = ' + repr(self.fwd.f))
            plt.savefig(self.outDir + 'Figs/fig76')
        