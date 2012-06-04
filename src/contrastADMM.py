'''
Created on Jun 3, 2012

@author: dstrauss

implementation of contrast source ADMM optimization
'''

import numpy as np
from maxwell import twoDim
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import sparseTools as spTools
import scipy.io as spio


class problem(twoDim):
    ''' class that extents the contrast - Xadmm algorithm '''
    def initOpt(self,rho,xi,uHat):
        self.rho = rho
        self.xi = xi
        self.uHat = uHat
        
        # add some local vars for ease
        self.s = 1j*self.muo*self.w
        self.A = self.nabla2+self.getk(0)
        
        # contrast variable
        self.X = np.zeros(self.nRx*self.nRy)
        # dual variable 
        self.Z = np.zeros(self.nRx*self.nRy)
        # scattered fields
        self.us = np.zeros(self.nx*self.ny)
        # just to make life easier:
        self.ub = self.sol[0].flatten()
        
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        self.Z = self.Z + (self.X - (self.s*self.Md*(self.ub + self.us))*P)
        
        uHatLocal =  self.uHat - self.Ms*self.ub  #remove background field
        
        pm = sparse.spdiags(self.s*P, 0, self.nRx*self.nRy, self.nRx*self.nRy)
        ds = pm*self.Md #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        bmuu = self.Ms.T*self.Ms + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmul = self.A.T.conj()
        
        rhsu = self.Ms.T.conj()*uHatLocal - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        
  
        bmxu = -self.rho*ds
        bmxx = self.rho*sparse.eye(self.nRx*self.nRy, self.nRx*self.nRy)
        bmxl = self.Md
        rhsx = self.rho*ds*self.ub - self.rho*self.Z
  
        bmlu = self.A
        bmlx = self.Md.T.conj()

        bmll = sparse.coo_matrix((self.nx*self.ny, self.nx*self.ny)) 
        rhsl = np.zeros(self.nx*self.ny)
  
  
        bm = spTools.vCat([spTools.hCat([bmuu, bmux, bmul]), \
                             spTools.hCat([bmxu, bmxx, bmxl]), \
                             spTools.hCat([bmlu, bmlx, bmll])])
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        
        updt = lin.spsolve(bm.tocsr(), rhsbm)
        
        N = self.nx*self.ny
        self.us = updt[:N]
        self.X = updt[N:(N+self.nRx*self.nRy)]
        
    def writeOut(self, ind=0):
        D = {'f':self.f, 'angle':self.incAng, 'sigMat':self.sigmap[0], 'ub':self.sol[0], \
             'us':self.us.reshape(self.nx,self.ny), 'uTrue':self.sol[1], \
             'X':self.X.reshape(self.nRx,self.nRy)}
        
        spio.savemat('contrastX' + repr(self.rank), D)
        
        
    
def aggregatorSerial(S, lmb, uBound):
    ''' routine to do the aggregation step and return an update for P '''
    N = np.size(S)
    n = S[0].nRx*S[0].nRy
    # print N
    # print n
    
    U = np.zeros((n,N),dtype='complex128')
    Q = np.zeros((n,N),dtype='complex128')
    
    for ix in range(N):
        s = S[ix].s
        U[:,ix] = s*S[ix].Md*(S[ix].ub + S[ix].us)
        Q[:,ix] = S[ix].X + S[ix].Z
        
    num = np.sum(U.real*Q.real + U.imag*Q.imag,1)
    
    den = np.sum(U.conj()*U,1) + lmb/S[0].rho
    
    P = (num/den).real
    P = np.maximum(P,0)
    P = np.minimum(P,uBound)
    
    gap = np.zeros(N)
    for ix in range(N):
        gap[ix] = np.linalg.norm(U[:,ix]*P - S[ix].X)
        
    return P

def aggregatorParallel(S, lmb, uBound, comm):
    ''' Do the aggregation step in parallel whoop! '''
    N = np.size(S)
    print repr(N) + ' == better be 1!'
    
    U = S.s*S.Md*(S.ub+S.us)
    Q = S.X + S.Z
    
    q = U.real*Q.real + U.imag*Q.imag
    num = np.zeros(q.shape)
    
    num = comm.allreduce(q,num,op=MPI.SUM)
    
    q = U.conj()*U + lmb/S.rho
    den = np.zeros(q.shape)
    den = comm.allreduce(q,den,op=MPI.SUM)
    
    P = (num/den).real
    P = np.maximum(P,0)
    P = np.minimum(P,uBound)
    
    gap = np.linalg.norm(U*P - S.X)
    print 'Proc ' + repr(comm.Get_rank()) + ' gap = ' + repr(gap)
    
    return P


def plotSerial(S,P,resid):
    import matplotlib.pyplot as plt
    ''' plotting routine for the serial passes '''
    N = np.size(S)
    plt.figure(383)
    plt.plot(resid)
    plt.savefig('fig383')
    
    plt.figure(387)
    plt.imshow(P.reshape(S[0].nRx,S[0].nRy), interpolation='nearest')
    plt.colorbar()
    plt.savefig('fig387')
    
    for ix in range(N):
        plt.figure(50+ix)
        # vv = S[ix].Ms*S[0].v
        uu = S[ix].Ms*S[0].us
        ub = S[ix].Ms*S[0].ub
        skt = S[ix].uHat-ub
    
        # io.savemat('uHat'+repr(ix), {'uh':uHat, 'ub':ub, 'skt':skt})

        # plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real, np.arange(S[0].nSen), vv.real)
        plt.plot(np.arange(S[0].nSen), skt.real, np.arange(S[0].nSen), uu.real)
        plt.savefig('fig'+repr(50+ix))
        
    plt.figure(76)
    plt.subplot(121)
    plt.imshow(S[0].us.reshape(S[0].nx,S[0].ny).real)
    plt.colorbar()
    
    plt.subplot(122)
    plt.imshow(S[1].us.reshape(S[0].nx,S[0].ny).real)
    plt.colorbar()
    plt.savefig('fig76')
    # plt.show()
    
def plotParallel(S,P,resid,rank):
    import matplotlib
    matplotlib.use('PDF')
    import matplotlib.pyplot as plt
    ''' Plotting routine if things are parallel'''
    # vv = S.Ms*S.v
    uu = S.Ms*S.us
    ub = S.Ms*S.ub
    skt = S.uHat-ub
    
    plt.figure(100+rank)
    plt.plot(np.arange(S.nSen), skt.real, np.arange(S.nSen), uu.real)
    plt.savefig('fig' + repr(100+rank))
    
    if rank==0:
        # then print some figures   
        plt.figure(383)
        plt.plot(resid)
        plt.savefig('fig383')
    
        plt.figure(387)
        plt.imshow(P.reshape(S.nRx,S.nRy), interpolation='nearest')
        plt.colorbar()
        plt.savefig('fig387')

        plt.figure(76)
        plt.subplot(121)
        plt.imshow(S.us.reshape(S.nx,S.ny).real)
        plt.colorbar()
    
        plt.subplot(122)
        plt.imshow(S.us.reshape(S.nx,S.ny).imag)
        plt.colorbar()
        plt.savefig('fig76')
    
    # all show!
    # plt.show()
