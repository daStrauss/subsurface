'''
Created on Jun 4, 2012

@author: dstrauss
'''

import maxwell.twoDim as twoDim
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
import sparseTools as spt
import numpy as np

class problem(twoDim):
    '''a class to do the born approximation iterations '''
    
    #internal parameters of the solver
    eabs = 1e-4
    erel = 1e-3
    rho = 0.0001
    maxit = 2000
    
    def initOpt(self,uHat, stepSize=0.005, upperBound=0.05):
        self.stepSize = stepSize
        self.uHat = uHat
        self.s = self.w*self.muo*1j
        self.uppB = upperBound
        
        
    def trueSolve(self,P):
        ''' an "internal" method to use for obtaining the current value of us '''
        A = self.nabla2 + self.getk(0) + sparse.spdiags(self.s*self.Md.T*P,0,self.nx*self.ny, self.nx*self.ny)
        
        return lin.spsolve(A,self.rhs,use_Umfpack=True),A
    
    def runOpt(self,P):
        ''' runtime module'''
        # get the local background -- requires a new solve!
        ub,A = self.trueSolve(P)
        localuHat = self.uHat - self.Ms*ub
        
        # produce some local matrices
        B = self.s*sparse.spdiags(ub,0,self.nx*self.ny,self.nx*self.ny)*self.Md.T
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
        
        M = spt.vCat([spt.hCat([Z.T*Z, sparse.coo_matrix((n,m)), A.T.conj()]), \
                      spt.hCat([sparse.coo_matrix((m,n)), self.rho*sparse.eye(m,m), B.T.conj()]),\
                      spt.hCat([A, B, sparse.coo_matrix((n,n))])])
        
        Q = lin.factorized(M)
        
        lowb = -P.copy()
        lowb[P-self.stepSize > 0 ] = -self.stepSize
        
        uppb = self.uppB-P
        uppb[P+self.stepSize < self.uppB] = self.stepSize
        
        okgo = True
        iterk = 0
        
        while okgo:
            iterk += 1
            rhs = np.concatenate((Z.T*localuHat, self.rho*(q-r), c))
            updt = Q(rhs)
            
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
    
    
    def aggregatorSerial(self,S,alpha):
        ''' super simple aggregator for sba method'''
        N = np.size(S)
        n = S[0].nRx*S[0].nRy
        
        P = np.zeros(n)
        
        for ix in range(N):
            P += S[ix].deltaP
        
        return (1.0/N)*P*alpha
    
    def aggregatorParallel(self, alpha,comm):
        pass
    
    def plotParallel(self):
        pass
    
    def plotSerial(self,S):
        pass
        
            