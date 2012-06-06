'''
Created on Jun 4, 2012

@author: dstrauss
'''

from maxwell import twoDim
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
import sparseTools as spt
import numpy as np
from mpi4py import MPI
import scipy.io as spio
import gc

class problem(twoDim):
    '''a class to do the born approximation iterations '''
    
    #internal parameters of the solver
    eabs = 1e-4
    erel = 1e-3
    rho = 0.0001
    maxit = 2000
    
    def initOpt(self, uHat, rho=0.005, xi=0.9, uBound=0.05, lmb=0):
        self.stepSize = rho
        self.alpha = xi
        self.uHat = uHat
        self.s = self.w*self.muo*1j
        self.uBound = uBound
        self.lmb = lmb
        
        
    def trueSolve(self,P):
        ''' an "internal" method to use for obtaining the current value of us '''
        
        self.A = self.nabla2 + self.getk(0) + sparse.spdiags(self.s*self.Md.T*P,0,self.nx*self.ny, self.nx*self.ny)
        self.us = lin.spsolve(self.A,self.rhs)

    
    def runOpt(self,P):
        ''' runtime module'''
        # get the local background -- requires a new solve!
        self.trueSolve(P)
        localuHat = self.uHat - self.Ms*self.us
        
        # produce some local matrices
        B = self.s*sparse.spdiags(self.us,0,self.nx*self.ny,self.nx*self.ny)*self.Md.T
        c = np.zeros(self.nx*self.ny)
        Z = self.Ms
        
        # grab some internal dimensions -- non ROM for now.
        n = self.nx*self.ny
        m = self.nRx*self.nRy
        
        # initialize some empty variables
        v = np.zeros(n) # update to fields (du)
        p = np.zeros(m) # estimate 1 for materials
        q = np.zeros(m) # estimate 2 for materials
        r = np.zeros(m) # dual variable for materials
        
        M = spt.vCat([spt.hCat([Z.T*Z, sparse.coo_matrix((n,m)), self.A.T.conj()]), \
                      spt.hCat([sparse.coo_matrix((m,n)), self.rho*sparse.eye(m,m), B.T.conj()]),\
                      spt.hCat([self.A, B, sparse.coo_matrix((n,n))])]).tocsc()
        
        print 'M format ' + repr(M.format)
        Q = lin.splu(M)
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
            updt = Q.solve(rhs)
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
                print 'At ADMM internal limit at iter ' + repr(iterk) 
                okgo = False
            elif (iterk >= self.maxit):
                okgo = False
                print 'Hit MAXX iterations'
            else:
                okgo = True
                
            r = r + (p-q)
        
        
        self.deltaP = q
        return q
    
    
    def aggregatorSerial(self,S):
        ''' super simple aggregator for sba method'''
        N = np.size(S)
        n = S[0].nRx*S[0].nRy
        
        P = np.zeros(n)
        
        for ix in range(N):
            P += S[ix].deltaP
        
        return (1.0/N)*P*self.alpha
    
    def aggregatorParallel(self,comm):
        ''' SUPER simple aggregator for sba method for working over parallel '''
        dP = np.zeros(self.nRx*self.nRy)
        
        dP = comm.allreduce(self.deltaP,dP,op=MPI.SUM)
        dP = dP*(1.0/comm.Get_size())*self.alpha
        return dP
    
    def plotParallel(self,P,resid,rank):
        ''' Plotting routine if things are parallel'''
        import matplotlib
        matplotlib.use('PDF')
        import matplotlib.pyplot as plt
        import os
        
        if not os.path.exists('sbaFigs'):
            os.makedirs('sbaFigs')
            
        
        # vv = S.Ms*S.v
        self.trueSolve(P)
            
        uu = self.Ms*(self.us - self.sol[0].flatten())
        ub = self.Ms*(self.sol[0].flatten())
        skt = self.uHat-ub
        
        plt.figure(100+rank)
        plt.plot(np.arange(self.nSen), skt.real, np.arange(self.nSen), uu.real)
        plt.savefig('sbaFigs/fig' + repr(100+rank))
        
        if rank==0:
            # then print some figures   
            plt.figure(383)
            plt.plot(resid)
            plt.savefig('sbaFigs/fig383')
        
            plt.figure(387)
            plt.imshow(P.reshape(self.nRx,self.nRy), interpolation='nearest')
            plt.colorbar()
            plt.savefig('sbaFigs/fig387')
    
            plt.figure(76)
            plt.subplot(121)
            plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).real)
            plt.colorbar()
        
            plt.subplot(122)
            plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).imag)
            plt.colorbar()
            plt.title('Final Scattered Fields f = ' + repr(self.f))
            plt.savefig('sbaFigs/fig76')
        
    
    def plotSerial(self,S,P,resid):
        ''' plotting routine for the serial passes '''
        import matplotlib
        matplotlib.use('PDF')
        import matplotlib.pyplot as plt
        import os
        
        if not os.path.exists('sbaFigs'):
            os.makedirs('sbaFigs')
        
        N = np.size(S)
        plt.figure(383)
        plt.plot(resid)
        plt.savefig('sbaFigs/fig383')
        
        plt.figure(387)
        plt.imshow(P.reshape(S[0].nRx,S[0].nRy), interpolation='nearest')
        plt.colorbar()
        plt.savefig('sbaFigs/fig387')
        
        for ix in range(N):
            plt.figure(50+ix)
            # vv = S[ix].Ms*S[0].v
            S[ix].trueSolve(P)
            
            uu = S[ix].Ms*(self.us - S[ix].sol[0].flatten())
            ub = S[ix].Ms*(S[0].sol[0].flatten())
            skt = S[ix].uHat-ub
        
            # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})
    
            # plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real, np.arange(S[0].nSen), vv.real)
            plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real)
            plt.savefig('sbaFigs/fig'+repr(50+ix))
            
        plt.figure(76)
        plt.subplot(121)
        plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).real)
        plt.colorbar()
    
        plt.subplot(122)
        plt.imshow((self.us.reshape(self.nx,self.ny)-self.sol[0]).imag)
        plt.colorbar()
        plt.savefig('sbaFigs/fig76')
        plt.title('Final Scattered Fields f = ' + repr(self.f))

        plt.savefig('sbaFigs/fig76')
        # plt.show()
    def writeOut(self):
        '''routine to print out information about the solve '''
        import os
        if not os.path.exists('sbaData'):
            os.mkdir('sbaData')
        
        D = {'f':self.f, 'angle':self.incAng, 'sigMat':self.sigmap[0], 'ub':self.sol[0], \
             'us':self.us.reshape(self.nx,self.ny), 'uTrue':self.sol[1]}
        
        spio.savemat('sbaData/sba' + repr(self.rank), D)
    
            