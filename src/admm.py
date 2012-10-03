'''
Created on May 24, 2012

@author: dstrauss

This set of routines implements the fast field splitting technique. 
Both parallel and serial methods are implemented
'''

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import scipy.io as spio
from optimize import optimizer

class problem(optimizer):
    '''A function for implementing the field splitting methods'''
#    def __init__(self, freq, rho, xi):
#        self.f=freq
#        self.rho = rho
#        self.xi = xi
    
    def initOpt(self,uHat, D):
        ''' prepare for upcoming iterations '''
        self.rho = D['rho']
        self.xi = D['xi']
        self.uHat = uHat
        self.lmb = D['lmb']
        self.uBound = D['uBound']
        self.F = np.zeros(self.fwd.N,dtype='complex128')
        self.E = np.zeros(self.fwd.N,dtype='complex128')
        self.us = np.zeros(self.fwd.N,dtype='complex128')
        self.v = np.zeros(self.fwd.N,dtype='complex128')
        self.ub = self.fwd.sol[0] # -- update sol should be flat
        self.obj = np.zeros(D['maxIter'])
        
        # the update for the first step can be precomputed
        self.A = self.fwd.nabla2+self.fwd.getk(0)
        self.s = self.fwd.getS()
#        print self.xi
        self.Q = sparse.vstack([self.A, -self.xi*sparse.eye(self.fwd.N,self.fwd.N)])
        self.Moo = self.rho*(self.Q.conj().T*self.Q) + self.fwd.Ms.T*self.fwd.Ms
        
        self.M = lin.factorized(self.Moo.tocsc())

    def runOpt(self, P):
        ''' method to run in each interation of the ADMM routine - super parallel '''
        
        #update dual variables
        self.F = self.F + (self.A*self.us + self.s*(self.fwd.Md.T*P)*(self.ub+self.v))
        self.E = self.E + self.xi*(self.v-self.us)
        
        # reasses local uHat
        uHatLocal = self.uHat - (self.fwd.Ms*self.ub)
        
        #solve for us update
        b = self.fwd.Ms.T*uHatLocal - \
        self.rho*(self.Q.conj().T*np.append((self.s*(self.ub+self.v)*(self.fwd.Md.T*P) + self.F), (self.xi*self.v)+self.E))
        self.us = self.M(b)
        
        # I was curious about solver tolerance -- bottom line is that it is good
        # print np.linalg.norm(self.Moo*self.us - b)
        
        # solve for v update
        Q = sparse.vstack([sparse.spdiags((self.s*self.fwd.Md.T*P),0,self.fwd.N,self.fwd.N), self.xi*sparse.eye(self.fwd.N,self.fwd.N)])
        M = self.rho*(Q.conj().T*Q)
        b = -self.rho*(Q.conj().T*np.append(self.A*self.us + self.F + self.s*(self.fwd.Md.T*P)*self.ub, (self.E - self.xi*self.us)))
        
        self.v = lin.spsolve(M,b)
        
        #oh that's awkward. obj has to get passed out in order to be indexed.
        obj = np.linalg.norm(self.fwd.Ms*self.v-uHatLocal)
        return obj
        # gap = np.linalg.norm(self.v-self.us)
        # print 'obj = ' + repr(obj) + ' Split Gap ' + repr(gap)
        # return obj
    
    def writeOut(self, rank, ix=0):
        import os
        if not os.path.exists(self.outDir + 'Data'):
            os.mkdir(self.outDir + 'Data')
        #oh awkward again -- need to reshape and so on so that things get written out in the right column
        # form.
        sgm = self.fwd.parseFields(self.fwd.sigmap[0])
        ub = self.fwd.parseFields(self.fwd.sol[0])
        ut = self.fwd.parseFields(self.fwd.sol[1])
        
        us = self.fwd.parseFields(self.us)
        v = self.fwd.parseFields(self.v)
        
        D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgm[0], 'ub':ub[0], \
             'us':us[0], 'uTrue':ut[0], \
             'v':v[0], 'rho':self.rho, 'xi':self.xi, 'obj':self.obj,\
             'flavor':self.flavor}
    
        spio.savemat(self.outDir + 'Data/splitField' + repr(rank) + '_' + repr(ix), D)
        
  
    
    def aggregatorSerial(self, S):
        '''Do the aggregate step updates '''
        # currently assumes that Md*Md' = I
        N = np.size(S)
        n = S[0].fwd.nRx*S[0].fwd.nRy
        # print N
        # print n
        
        U = np.zeros((n,N),dtype='complex128')
        Q = np.zeros((n,N),dtype='complex128')
        
        for ix in range(N):
            s = S[ix].s
            U[:,ix] = s*S[ix].fwd.Md*(S[ix].ub+S[ix].v)
            Q[:,ix] = S[ix].fwd.Md*((S[ix].A*S[ix].us) + S[ix].F)
            
        num = np.sum(U.real*Q.real + U.imag*Q.imag,1)
        print num.shape
        den = np.sum(U.conj()*U,1) + self.lmb/S[0].rho
        
        P = (-num/den).real
        P = np.maximum(P,0)
        P = np.minimum(P,self.uBound)
        
        gap = np.zeros(N)
        for ix in range(N):
            gap[ix] = np.linalg.norm(S[ix].fwd.Md.T*(U[:,ix]*P) + S[ix].A*S[ix].us)
            
        return P

    def aggregatorParallel(self,comm):
        ''' Do the aggregation step in parallel whoop! '''
        # currently assumes that Md*Md' is I
        N = np.size(self)
        print repr(N) + ' == better be 1!'
        
        U = self.s*self.fwd.Md*(self.ub+self.v)
        Q = self.fwd.Md*((self.A*self.us) + self.F)
        
        q = U.real*Q.real + U.imag*Q.imag
        num = np.zeros(q.shape)
        
        num = comm.allreduce(q,num,op=MPI.SUM)
        
        q = U.conj()*U + self.lmb/self.rho
        den = np.zeros(q.shape)
        den = comm.allreduce(q,den,op=MPI.SUM)
        
        P = (-num/den).real
        P = np.maximum(P,0)
        P = np.minimum(P,self.uBound)
        
        gap = np.linalg.norm(self.fwd.Md.T*(U*P) + self.A*self.us)
        print 'Proc ' + repr(comm.Get_rank()) + ' gap = ' + repr(gap)
        
        return P

    def aggregatorSemiParallel(self,S,comm):
        ''' Do the aggregation step in parallel whoop! '''
        N = np.size(S)
        n = S[0].fwd.nRx*S[0].fwd.nRy
        
        QL = sparse.lil_matrix( (n,n), dtype='complex128' )
        bL = np.zeros(n, dtype='complex128')
        # U = np.zeros((n,N),dtype='complex128')
        # Q = np.zeros((n,N),dtype='complex128')
        
        for L in S:
            D = sparse.spdiags(L.s*(L.ub+L.v),0, self.fwd.N,self.fwd.N)
            g = D*L.fwd.Md.T
            bL += g.T.conj()*((L.A*L.us) + L.F)
            QL += g.T.conj()*g

#        QL = QL/N
#        bL = bL/N
        
        Q = sparse.lil_matrix((n,n), dtype='complex128')
        
        Q = comm.allreduce(QL,Q,op=MPI.SUM)
#        Q = Q/comm.Get_size()
        
        b = np.zeros(n)
        b = comm.allreduce(bL,b,op=MPI.SUM)
#        b = b/comm.Get_size()
        
        # hah, no 1/n's because the all got em
        Q = Q + self.lmb*sparse.eye(n,n)
        P = lin.spsolve(Q,-b).real
        P = np.maximum(P,0)
        P = np.minimum(P,self.uBound)
        
        # gap = np.linalg.norm(self.Md.T*(U*P) + self.A*self.us)
        # print 'Proc ' + repr(comm.Get_rank()) + ' gap = ' + repr(gap)
        
        return P

    def plotSerial(self,S,P,resid):
        ''' plotting routine for the serial passes '''
        import matplotlib
        matplotlib.use('PDF')
        import matplotlib.pyplot as plt
        plt.close('all')
        import os
        if not os.path.exists('Figs'):
            os.mkdir('Figs')
        
        N = np.size(S)
        plt.figure(383)
        plt.plot(resid)
        plt.savefig(self.outDir + 'Figs/fig383')
        
        plt.figure(387)
        plt.imshow(P.reshape(S[0].nRx,S[0].nRy), interpolation='nearest')
        plt.colorbar()
        plt.savefig(self.outDir + 'Figs/fig387')
        
        for ix in range(N):
            plt.figure(50+ix)
            vv = S[ix].fwd.Ms*S[0].v
            uu = S[ix].fwd.Ms*S[0].us
            ub = S[ix].fwd.Ms*S[0].ub
            skt = S[ix].uHat-ub
        
            # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})
    
            plt.plot(np.arange(S[0].fwd.nSen), skt.real, np.arange(S[0].fwd.nSen), uu.real, np.arange(S[0].fwd.nSen), vv.real)
            plt.savefig(self.outDir + 'Figs/fig'+repr(50+ix))
        
        u = self.fwd.parseFields(S[0].us)    
        plt.figure(76)
        plt.subplot(121)
        plt.imshow(u[0].real)
        plt.colorbar()
        
        lv = self.fwd.parseFields(S[0].v)
        
        plt.subplot(122)
        plt.imshow(lv.real)
        plt.colorbar()
        # plt.show()
        plt.savefig(self.outDir + 'Figs/fig76')
    
    def plotParallel(self,P,resid,rank):
        ''' Plotting routine if things are parallel'''
        import matplotlib
        matplotlib.use('PDF')
        import matplotlib.pyplot as plt
        plt.close('all')
        import os
        
        if not os.path.exists('Figs'):
            os.mkdir('Figs')
        
        vv = self.fwd.Ms*self.v
        uu = self.fwd.Ms*self.us
        ub = self.fwd.Ms*self.ub
        skt = self.uHat-ub
        
        plt.figure(100+rank)
        plt.plot(np.arange(self.fwd.nSen), skt.real, np.arange(self.fwd.nSen), uu.real, np.arange(self.fwd.nSen), vv.real)
        plt.savefig(self.outDir + 'Figs/fig' + repr(100+rank))
        
        if rank==0:
            # then print some figures   
            plt.figure(383)
            plt.plot(resid)
            plt.savefig(self.outDir + 'Figs/fig383' )
        
            plt.figure(387)
            plt.imshow(P.reshape(self.fwd.nRx, self.fwd.nRy), interpolation='nearest')
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig387'  )
    
            u = self.fwd.parseFields(self.us)
            plt.figure(76)
            plt.subplot(121)
            plt.imshow(u[0].real)
            plt.colorbar()
            
            vl = self.fwd.parseFields(self.v)
            plt.subplot(122)
            plt.imshow(vl[0].real)
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig76'  )
        
        # all show!
#        plt.show()
        
 
    
    def plotSemiParallel(self,P,resid,rank,ix=0):
        ''' Plotting routine if things are parallel'''
        import os
        import matplotlib.pyplot as plt
        plt.close('all')
        
        if not os.path.exists(self.outDir+'Figs'):
            os.mkdir(self.outDir + 'Figs')
        
        vv = self.fwd.Ms*self.v
        uu = self.fwd.Ms*self.us
        ub = self.fwd.Ms*self.ub
        skt = self.uHat-ub
        
        plt.figure(100+rank+10*ix)
        plt.plot(np.arange(self.fwd.nSen), skt.real, np.arange(self.fwd.nSen), uu.real, np.arange(self.fwd.nSen), vv.real)
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
            
            vl = self.fwd.parseFields(self.v)
            plt.subplot(122)
            plt.imshow(vl[0].real)
            plt.colorbar()
            plt.savefig(self.outDir + 'Figs/fig76'  )
        
        # all show!
#        plt.show()
        

    
    
    
    
    
    