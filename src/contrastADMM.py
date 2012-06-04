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


class problem(twoDim):
    ''' class that extents the contrast - Xadmm algorithm '''
    def initOpt(self,rho,xi,uHat):
        self.rho = rho
        self.xi = xi
        self.uHat = uHat
        self.s = 1j*self.muo*self.w
        
        # contrast variable
        self.X = np.zeros(self.nRx*self.nRy)
        # dual variable 
        self.Z = np.zeros(self.nRx*self.nRy)
        # scattered fields
        self.us = np.zeros(self.nx*self.ny)
        # just to make life easier:
        self.u = self.sol[0].flatten()
        
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        self.Z = self.Z + (self.X - (self.s*self.Md*(self.u + self.us))*P)
        
        uHatLocal =  self.uHat - self.Ms*self.u  #remove background field
        
        pm = sparse.spdiags(self.s*P, 0, self.nRx*self.nRy, self.nRx*self.nRy)
        ds = pm*self.Md #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        bmuu = self.Ms.T*self.Ms + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmul = self.A.T.conj()
        
        rhsu = self.Ms.T.conj()*uHatLocal - self.rho*(ds.T.conj()*ds)*self.u + self.rho*ds.T.conj()*self.Z
        
  
        bmxu = -self.rho*ds
        bmxx = self.rho*sparse.eye(self.nRx*self.nRy, self.nRx*self.nRy)
        bmxl = self.Md
        rhsx = self.rho*ds*self.u - self.rho*self.Z
  
        bmlu = self.A
        bmlx = F.Md.T.conj()
        bmll = np.empty(F.nx*F.nx,F.nx*F.nx,0);
        
        rhsl = np.zeros(self.nx*self.ny)
  
  
        bm = spTools.spVcat([spTools.spHcat([bmuu, bmux, bmul]), \
                             spTools.spHcat([bmxu, bmxx, bmxl]), \
                             spTools.spHcat([bmlu, bmlx, bmll])])
        
        rhsbm = concatenate((rhsu, rhsx, rhsl))
        
  
    updt = bm\rhsbm;
    F.opt.us = updt(1:F.nx*F.nx);
    F.opt.X = updt((F.nx*F.nx+1):(F.nx*F.nx+F.nRx*F.nRy));
    
def aggregatorSerial(S, P, lmb, uBound):
    pass

def aggregatorParllel(S, P, lmb, uBound, comm):
    pass
