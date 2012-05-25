'''
Created on May 24, 2012

@author: dstrauss
'''

import numpy as np
from maxwell import twoDim
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
import numpy

class fieldSplit(twoDim):
    '''A function for implementing the field splitting methods'''
#    def __init__(self, freq, rho, xi):
#        self.f=freq
#        self.rho = rho
#        self.xi = xi
    
    def initOpt(self,rho,xi,uHat):
        ''' prepare for upcoming iterations '''
        self.rho = rho
        self.xi = xi
        self.uHat = uHat
        self.F = np.zeros(self.N,dtype='complex128')
        self.E = np.zeros(self.N,dtype='complex128')
        self.us = np.zeros(self.N,dtype='complex128')
        self.v = np.zeros(self.N,dtype='complex128')
        self.ub = self.sol[0].flatten()
        
        # the update for the first step can be precomputed
        self.A = self.nabla2+self.getk(0)
        
        self.Q = sparse.vstack([self.A, -self.xi*sparse.eye(self.N,self.N)])
        M = self.rho*(self.Q.T*self.Q) + self.Ms.T*self.Ms
        
        self.M = lin.factorized(M.tocsc())

    def runOpt(self, P):
        ''' method to run in each interation of the ADMM routine - super parallel '''
        
        #update dual variables
        self.F = self.F + (self.A*self.us + self.c*(self.Md.T*P)*(self.ub+self.v))
        self.E = self.E + self.xi*(self.v-self.us)
        
        # reasses local uHat
        uHatLocal = self.uHat - (self.Ms*self.ub)
        
        #solve for us update
        b = self.Ms.T*uHatLocal - self.rho*(self.Q.T*np.append(self.c*(self.ub+self.v)*(self.Md.T*P), self.xi*self.v+self.E))
        self.us = self.M(b)
        
        # solve for v update
        Q = sparse.vstack([sparse.spdiags((self.c*self.Md.T*P),0,self.N,self.N), self.xi*sparse.eye(self.N,self.N)])
        M = self.rho*(Q.T*Q)
        b = -self.rho*(Q.T*np.append(self.A*self.us + self.F + self.c*(self.Md.T*P), self.E - self.xi*self.us))
        
        self.v = lin.spsolve(M,b)
        
        obj = np.linalg.norm(self.Ms*self.v-uHatLocal)
        gap = np.linalg.norm(self.v-self.us)
        return obj, gap
    
    
    
def aggregateFS(S,lmb,uBound):
    '''Do the aggregate step updates '''
    N = np.size(S)
    n = S[0].nRx*S[0].nRy
    
    # U = np.zeros((n,N),dtype='complex128')
    # Q = np.zeros((n,N),dtype='complex128')
    
    for ix in range(N):
        c = S[ix].c
        U = np.hstack([U,(S[ix].Md*S[ix].ub+S[ix].Md*S[ix].v)])
        Q[:,ix] = S[ix].Md*((S[ix].A*S[ix].us) + S[ix].F)
        
    num = np.sum(U.real*Q.real + U.imag*Q.imag,1)
    den = np.sum(U.conj()*U,1) + lmb/S[0].rho
    
    P = (-num/den).real
    P = max(P,0)
    P = min(P,uBound)
    
    gap = np.zeros(N)
    for ix in range(N):
        gap[ix] = np.linalg.norm(S[ix].Md.T*(U[:,ix]*P) + S[ix].A*S[ix].us)
        
    return P,gap
    
    