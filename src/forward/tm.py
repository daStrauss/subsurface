'''
Created on Jun 12, 2012

@author: dstrauss
'''

from model import fwd
from scipy import sparse
import sparseTools as spt
import numpy as np

class solver(fwd):
    def setspace(self, nx,ny,dx,dy):
        super(solver,self).setspace(nx,ny,dx,dy)
        # self.ExNx = nx+1
        # self.ExNy = ny
        # self.EyNx = nx
        # self.EyNy = ny+1
        
    def makeGradOper(self):
        hz2ex = sparse.kron(self.po*self.d2,sparse.eye(self.nx+1,self.nx+1));
        hz2ey = sparse.kron(sparse.eye(self.nx+1,self.nx+1), -self.po*self.d2);

        ex2hz = sparse.kron(-self.ph*self.d1,sparse.eye(self.nx+1,self.nx+1));
        ey2hz = sparse.kron(sparse.eye(self.nx+1,self.nx+1),self.ph*self.d1);

        n = self.nx*(self.nx+1);
        N = (self.nx+1)**2;
        
        self.nabla2 = spt.vCat([spt.hCat([sparse.coo_matrix((n,n)), sparse.coo_matrix((n,n)), hz2ex]), \
                                spt.hCat([sparse.coo_matrix((n,n)), sparse.coo_matrix((n,n)), hz2ey]), \
                                spt.hCat([ex2hz, ey2hz, sparse.coo_matrix((N,N))]) ])
                              
    def setMs(self, nSensors=10):
        '''Tell me the number of sensors, and I will distribute them equally across the surface
        '''
        self.nSen = nSensors
        indx = np.round(np.linspace(self.npml+10,self.nx-self.npml-10, nSensors)-1).astype(int);
        oprx = np.zeros((self.nx,self.ny),dtype='bool')
        
        oprx[indx,self.div] = 1;
        
        idx = np.arange(self.N)
        oprx = oprx.flatten()
        idx = idx[oprx]
        
        self.Ms = sparse.dok_matrix((self.N,idx.size))
        
        for i in range(sum(oprx)):
            self.Ms[idx[i],i] = 1.0
        self.Ms = self.Ms.tocsc()
        self.Ms = self.Ms.T
        
                              
                              
        