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
        In the TM script, this is going to take Magnetic Fields only.
        '''
        self.nSen = nSensors
        indx = np.round(np.linspace(self.npml+10,self.nx-self.npml-10, nSensors)-1).astype(int);
        
        # in TM case -- going for Hz
        oprx = np.zeros((self.nx+1,self.ny+1),dtype='bool')
        
        oprx[indx,self.div] = 1;
        
        # stub to get to the actual indices, because we can't do boolean slices
        idx = np.arange((self.nx+1)*(self.nx+1))
        oprx = oprx.flatten()
        idx = idx[oprx]
        
        self.Ms = sparse.dok_matrix((self.N,idx.size))
        
        # there's probably a more pythonic way here ---
        for i in range(sum(oprx)):
            self.Ms[idx[i],i] = 1.0
        self.Ms = self.Ms.tocsc()
        self.Ms = self.Ms.T
        
        z = self.Ms.shape
        #super dimension of the ex,ey
        n = (self.nx+1)*self.nx
        self.Ms = spt.hCat([sparse.coo_matrix((z[0],n)), sparse.coo_matrix((z[0],n)), self.Ms])
        
    def setMd(self, xrng, yrng):
        '''Tell me the xrange and the yrange and I'll make selector
        '''
        # because we need to know -- only once
        self.nRx = xrng[1]-xrng[0]
        self.nRy = yrng[1]-yrng[0]
        
        oprx = np.zeros((self.nx+1,self.ny),dtype='bool')
        #might have to revise this later.
        oprx[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1 
        
        idx = np.arange((self.nx+1)*self.nx)
        oprx = oprx.flatten()
        idx = idx[oprx]
        Mdx = sparse.dok_matrix(((self.nx+1)*self.nx,idx.size), 'bool')
        
        for i in range(idx.size):
            Mdx[idx[i],i]=1
            
        opry = np.zeros((self.nx,self.ny+1),dtype='bool')
        # might have to revise this later.
        opry[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1
        
        idx = np.arange((self.nx+1)*self.nx)
        opry = opry.flatten()
        idx = idx[opry]
        Mdy = sparse.dok_matrix(((self.nx+1)*self.nx,idx.size), 'bool')
        
        for i in range(idx.size):
            Mdy[idx[i],i] = 1
        # only do one conversion to csc    
        N = (self.nx+1)*(self.nx+1)
        self.Md = spt.vCat([Mdx, Mdy, sparse.coo_matrix((N))]) #oh bugger, problems. 
        self.Md = self.Md.tocsc()
        self.Md = self.Md.T
                              
                              
        