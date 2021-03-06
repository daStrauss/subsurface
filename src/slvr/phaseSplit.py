'''
Created on Nov 11, 2012
Copyright © 2013
The Board of Trustees of The Leland Stanford Junior University.
All Rights Reserved

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@author: dstrauss


implementation of contrast source ADMM optimization
with phase split modification.

This is a new experiment to see if I can bring any of the power of the phase split method to
light here? 
So we know that X and u ought to have the same phase, but in the end, they rarely do. (odd? dontcha think?)
So lets add in the phase holding term. Cuz they ought to have the same phase.

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
        
        self.pp = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.X = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.Z = np.zeros(self.fwd.getXSize(),dtype='complex128')
        
        self.tD = np.zeros(self.fwd.nRx*self.fwd.nRy,dtype='complex128')
        self.tL = np.zeros(self.fwd.nRx*self.fwd.nRy,dtype='complex128')
        
        
        self.fwd.setCTRX()
        
        ''' subtract out the background field '''
        self.uHat = self.uHat - self.fwd.Ms*self.ub
        ''' in this instance, I don't care about the results, i.e., I don't care about the actual solutions'''
        
    def internalHard(self, thk):
        '''creates the matrix every time, a hard alternative, to internalSymbolic
        so that I can do A/B testing easily w.r.t the old standard'''
        nX = self.fwd.getXSize()
        pm = sparse.spdiags(self.s*self.fwd.p2x*thk, 0, nX, nX)
        plp = sparse.spdiags(np.exp(1j*self.pp),0,nX,nX)
        # print pm.shape
        # print self.fwd.x2u.shape
        ''' ds changes at ever iteration '''
        ds = pm*self.fwd.x2u.T #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        ''' changes '''
        bmuu = self.scaleC*(self.fwd.Ms.T*self.fwd.Ms) + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*(ds.T.conj()*plp)
        bmxu = -self.rho*(plp.T.conj()*ds)
        ''' static ''' 
        bmul = self.A.T.conj()
        bmxx = self.rho*sparse.eye(nX, nX)
        bmxl = plp.T.conj()*self.fwd.x2u.T
        bmlu = self.A
        bmlx = self.fwd.x2u*plp
        bmll = sparse.coo_matrix((self.fwd.N, self.fwd.N))
        
        ''' right hand side ''' 
        rhsu = self.fwd.Ms.T.conj()*self.uHat - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        rhsx = self.rho*plp.T.conj()*(ds*self.ub - self.Z) # chng
        rhsl = np.zeros(self.fwd.N)
  
  
        bm = spTools.vCat([spTools.hCat([bmuu, bmux, bmul]), \
                           spTools.hCat([bmxu, bmxx, bmxl]), \
                           spTools.hCat([bmlu, bmlx, bmll])])
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        
        updt = lin.spsolve(bm.tocsr(), rhsbm)
        
        # N = self.nx*self.ny
        us = updt[:self.fwd.N]
        x = updt[self.fwd.N:(self.fwd.N+nX)]
        return us,x
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        
        ''' the global tT update happens effectively here -- in parallel with u,x update'''
        
        ''' jointly update u,x '''
        # pfake = (self.upperBound/2.0)*np.ones(self.fwd.getXSize(),dtype='complex128')
        self.pp = np.angle(self.fwd.x2u.T*(self.us + self.ub))
        nX = self.fwd.getXSize()
        plp = sparse.spdiags(np.exp(1j*self.pp),0,nX,nX)
        
        self.us,self.X = self.internalHard(self.tL)
        
        ''' do some accounting '''
        self.gap.append(np.linalg.norm(plp*self.X - self.tL*(self.s*self.fwd.Md*(self.us+self.ub))))
        
        ''' update tL '''
        self.tL = self.updateThetaLocal(P)
        self.pL.append(self.tL)
        
        '''update dual variables last '''
        self.Z = self.Z + (plp*self.X - (self.s*self.fwd.x2u.T*(self.ub + self.us))*(self.fwd.p2x*self.tL))
        self.tD = self.tD + (self.tL-P)
        
        obj = np.linalg.norm(self.uHat-self.fwd.Ms*self.us)
        self.objInt.append(obj)
        return obj
    
    def updateThetaLocal(self,P):
        '''routine that solves for a local copy of the parameters theta '''
        uL = self.s*self.fwd.Md*(self.us + self.ub)
        
        nX = self.fwd.getXSize()
        plp = sparse.spdiags(np.exp(1j*self.pp),0,nX,nX)
        
        bL = plp*self.X + self.Z
            
        localT = (P - self.tD)
        
        theta = (self.rho*(uL.conj()*bL) + self.xi*(localT))/(self.rho*(uL.conj()*uL) + self.xi + self.lmb)
        
        theta = theta.real
        
        theta = np.maximum(theta,0)
        theta = np.minimum(theta,self.upperBound)
        return theta
        
    
    def writeOut(self, rank, ix=0):
        import os
        assert os.path.exists(self.outDir + 'Data')
        
        us = self.fwd.parseFields(self.us)
        ub = self.fwd.parseFields(self.fwd.sol[0])
        sgmm = self.fwd.parseFields(self.fwd.sigmap[0])
        uTrue = self.fwd.parseFields(self.fwd.sol[1])
            
        D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgmm[0], 'ub':ub[0], \
             'us':us[0], 'uTrue':uTrue[0], \
             'X':self.X, 'obj':self.obj, 'flavor':self.fwd.flavor, 'gap':self.gap, \
             'obj':self.objInt, 'Ms':self.fwd.Ms, 'phist':self.pL}
        
        spio.savemat(self.outDir + 'Data/contrastX' + repr(rank) + '_' + repr(ix), D)
        

    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! 
        Revised to be just a simple aggregation/mean step '''
        
        n = self.fwd.nRx*self.fwd.nRy
        tt = self.tL + self.tD
        
        T = np.zeros(n,dtype='complex128')
        
                
        T = comm.allreduce(tt,T,op=MPI.SUM)
        T = T/comm.Get_size()
                       
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
