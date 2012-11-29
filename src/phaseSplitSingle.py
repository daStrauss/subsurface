'''
Created on Nov 20, 2012

@author: dstrauss

copy out from phaseSplit:
I want to implement a single frequency, efficient method for the phase split method
Not going to worry about consensus for the moment. Maybe it works for 3d? (unlikely, but we'll try).
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

        ''' create some lists so that I can see how the algorithm progresses'''        
        self.gap = list()
        self.objInt = list()
        self.pL = list()
        self.phaseList = list()
        self.rlist = list()
        
        self.us = np.zeros(self.fwd.N,dtype='complex128')
        # just to make life easier:
        self.ub = self.fwd.sol[0] # shouldn't need --> .flatten()
        
        self.pp = np.zeros(self.fwd.getXSize(),dtype='complex128')
        
        self.r = np.zeros(self.fwd.getXSize(),dtype='complex128') # primal
        self.rt = np.zeros(self.fwd.getXSize(),dtype='complex128') # twiddle
        self.rd = np.zeros(self.fwd.getXSize(), dtype='complex128') # dual
        
        self.z = np.zeros(2*self.fwd.getXSize(), dtype='complex128') #primal?
        self.zt = np.zeros(self.fwd.getXSize(), dtype='complex128') # twiddle
        self.zd = np.zeros(self.fwd.getXSize(), dtype='complex128') # dual
        
        # self.pp = sparse.spdiags(np.exp(1j*np.angle(self.fwd.x2u.T*self.ub)),0,self.fwd.getXSize(),self.fwd.getXSize())
        
        self.fwd.setCTRX()
        
        ''' subtract out the background field '''
        self.uHat = self.uHat - self.fwd.Ms*self.ub
        
        ''' create the system KKT matrix '''
        uu = self.fwd.Ms.T*self.fwd.Ms + self.rho*(self.upperBound*self.fwd.x2u*self.fwd.x2u.T*self.upperBound)
        ur = -self.rho*self.upperBound*self.fwd.x2u
        ul = self.A.T.conj()
        
        rr = 2*self.rho*sparse.eye(self.fwd.getXSize(),self.fwd.getXSize())
        rl = self.s.conj()*self.fwd.x2u.T
        
        ll = sparse.coo_matrix((self.fwd.N,self.fwd.N))
        
        SM = spTools.vCat([spTools.hCat([uu,ur,ul]),\
                           spTools.hCat([ur.T.conj(), rr, rl]), \
                           spTools.hCat([ul.T.conj(), rl.T.conj(), ll])])
                           
        self.internalGo = lin.factorized(SM)
    
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        ''' jointly update u,r '''
        
        self.pp = np.angle(self.fwd.x2u.T*(self.us + self.ub))
        self.phaseList.append(self.pp)
        self.pL.append(self.fwd.x2u.T*self.us)
        
        nX = self.fwd.getXSize()
        plp = sparse.spdiags(np.exp(1j*self.pp),0,nX,nX)
        
        ru = self.fwd.Ms.T*self.uHat - self.rho*(self.upperBound*self.fwd.x2u*plp*(self.zt-self.zd)) - \
        self.rho*(self.upperBound*self.fwd.x2u*self.fwd.x2u.T*self.upperBound)*self.ub 
        
        # + something with ub
        rr = self.rho*plp*(self.zt-self.zd) + self.rho*plp*(self.rt-self.rd) + \
            self.rho*self.upperBound*self.fwd.x2u.T*(self.ub)
        rl = np.zeros(self.fwd.N)
        
        rhh = np.concatenate((ru,rr,rl))
        
        sol = self.internalGo(rhh)
        
        self.us = sol[:self.fwd.N]
        self.r = plp.T.conj()*sol[self.fwd.N:(self.fwd.N+nX)]
        
        self.rlist.append(self.r)
        
        self.rt = (self.r+self.rd).real
        self.rt = np.maximum(self.rt,0.0)
        
        self.zt = (self.r-self.upperBound*plp.T.conj()*self.fwd.x2u.T*(self.us+self.ub) + self.zd).real
        self.zt = np.minimum(self.zt,0.0)
        
        ''' update dual variables '''
        self.rd = self.rd + self.r-self.rt
        self.zd = self.zd + self.r-self.upperBound*plp.conj()*self.fwd.x2u.T*(self.us+self.ub) - self.zt
        
        
        gap = np.linalg.norm(self.r-self.upperBound*plp.conj()*self.fwd.x2u.T*(self.us+self.ub) - self.zt)
        obj = np.linalg.norm(self.uHat-self.fwd.Ms*self.us)
        self.objInt.append(obj)
        self.gap.append(gap)
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
             'rt':self.rt, 'rd':self.rd, 'zd':self.zd, 'zt':self.zt, \
             'obj':self.obj, 'flavor':self.fwd.flavor, 'gap':self.gap, \
             'obj':self.objInt, 'Ms':self.fwd.Ms, 'Md':self.fwd.Md,\
             'phist':self.pL,'pp':self.pp, 'phl':self.phaseList, 'rlist':self.rlist}
        
        spio.savemat(self.outDir + 'Data/contrastX' + repr(rank) + '_' + repr(ix), D)
        

    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! 
        Revised to be just a simple aggregation/mean step '''
        
        ''' nothing to do here really'''
        print self.rt.shape
        P = self.rt/np.abs(self.fwd.Md*(self.us+self.ub))
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
