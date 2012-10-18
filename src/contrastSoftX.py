'''
Created on Oct 17, 2012

@author: dstrauss

implementation of contrast source optimization with soft update - i.e. requires one factorization.

'''

import numpy as np
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI
import sparseTools as spt
import scipy.io as spio
from optimize import optimizer
# import matplotlib.pyplot as plt


class problem(optimizer):
    ''' an instance of an optimizer that will solve using the soft contrast X algorithm 
    In this algorithm, I'm going to discard the "background field" its just annoying.'''
    def initOpt(self, uHat, D):
        self.rho = D['rho']
        self.xi = D['xi']
        
        self.upperBound = D['uBound']
        self.lmb = D['lmb']
        self.obj = np.zeros(D['maxIter'])
        
        # add some local vars for ease
        self.s = self.fwd.getS() #  1j*self.muo*self.w
        self.A = self.fwd.nabla2+self.fwd.getk(0)
        
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
        '''remove background field once and for all'''
        self.uHat = uHat - self.fwd.Ms*self.ub
        # create some new operators for doing what is necessary for the 
        # contrast X work
        self.X = np.zeros(self.fwd.getXSize(),dtype='complex128')
        self.Z = np.zeros(self.fwd.getXSize(),dtype='complex128')
        
        # this operator creates a matrix that maps to the X variables which are field specific
        self.fwd.setCTRX()
        # create an operator that makes a projector used for updating u,x jointly
        self.prepProjector()
        
        
    
    def runOpt(self,P):
        ''' to run at each layer at each iteration '''
        
        # update dual variable
        self.Z = self.Z + (self.X - (self.s*self.fwd.x2u.T*(self.us+self.ub))*(self.fwd.p2x*P))
        
        
        
        # uHatLocal =  self.uHat - self.fwd.Ms*self.ub  #remove background field
        
        nX = self.fwd.getXSize()
        pm = sparse.spdiags(self.s*self.fwd.p2x*P, 0, nX, nX)
        # print pm.shape
        # print self.fwd.x2u.shape
        ds = pm*self.fwd.x2u.T #  The sampling and material scaling.
  
        # Construct the KKT Matrix
        bmuu = self.fwd.Ms.T*self.fwd.Ms + self.rho*(ds.T.conj()*ds)
        bmux = -self.rho*ds.T.conj()
        bmul = self.A.T.conj()
        
        rhsu = self.fwd.Ms.T.conj()*self.uHat - self.rho*(ds.T.conj()*ds)*self.ub + self.rho*ds.T.conj()*self.Z
        # 
  
        bmxu = -self.rho*ds
        bmxx = self.rho*sparse.eye(nX, nX)
        bmxl = self.fwd.x2u.T
        rhsx = - self.rho*self.Z + self.rho*ds*self.ub 
        # 
        
        bmlu = self.A
        bmlx = self.fwd.x2u

        bmll = sparse.coo_matrix((self.fwd.N, self.fwd.N)) 
        rhsl = np.zeros(self.fwd.N)
        # rhsl = self.fwd.rhs
  
        bm = spt.vCat([spt.hCat([bmuu, bmux, bmul]), \
                       spt.hCat([bmxu, bmxx, bmxl]), \
                       spt.hCat([bmlu, bmlx, bmll])])
        
        rhsbm = np.concatenate((rhsu, rhsx, rhsl))
        
        updt = lin.spsolve(bm.tocsr(), rhsbm)
        
        # N = self.nx*self.ny
        dirUS = updt[:self.fwd.N]
        dirX = updt[self.fwd.N:(self.fwd.N+nX)]
        
        
#        
        fooUs,fooX = self.contrastProjector(P,dirUS,dirX)
        
        D = {'dirUs':dirUS.reshape(self.N,self.N), 'dirX':dirX.reshape(self.fwd.nRx,self.fwd.nRy), 'fooUs':fooUs.reshape(self.N,self.N), 'fooX':fooX.reshape(self.fwd.nRx,self.nRy)}
        spio.savemat('dircomp',D)
        
        self.us = fooUs;
        self.X = fooX;
#        print 'udiff ' + repr(np.linalg.norm(usCP-self.us)/np.linalg.norm(self.us))
#        print 'xdiff ' + repr(np.linalg.norm(xCP-self.X)/np.linalg.norm(self.X))

        
        obj = np.linalg.norm(self.uHat-self.fwd.Ms*self.us)
        return obj
        
    def writeOut(self, rank, ix=0):
        import os
        assert os.path.exists(self.outDir + 'Data')
        
        us = self.fwd.parseFields(self.us)
        ub = self.fwd.parseFields(self.fwd.sol[0])
        sgmm = self.fwd.parseFields(self.fwd.sigmap[0])
        uTrue = self.fwd.parseFields(self.fwd.sol[1]-self.fwd.sol[0])
        
        uSmp = self.fwd.Ms*self.us
        uSolSmp = self.fwd.Ms*(self.fwd.sol[1]-self.fwd.sol[0])
            
        D = {'f':self.fwd.f, 'angle':self.fwd.incAng, 'sigMat':sgmm[0], 'ub':ub[0], \
             'us':us[0], 'uTrue':uTrue[0], \
             'X':self.X, 'obj':self.obj, 'flavor':self.fwd.flavor, 'uSmp':uSmp, 'uSolSmp':uSolSmp}
        
        spio.savemat(self.outDir + 'Data/contrastSoftX' + repr(rank) + '_' + repr(ix), D)


    def prepProjector(self):
        ''' create a projection function, stored inside, that does the u,x projection step'''
        n = self.fwd.N
        m = self.fwd.getXSize()
        
        uu = self.fwd.Ms.T*self.fwd.Ms+self.xi*sparse.eye(n,n)
        ux = sparse.coo_matrix((n,m))
        ul = self.A.T.conj()
        
        xu = sparse.coo_matrix((m,n))
        xx = self.xi*sparse.eye(m,m)
        xl = self.fwd.x2u.T.conj()
        
        lu = self.A
        lx = self.fwd.x2u
        ll = sparse.coo_matrix((n,n))
        
        M = spt.vCat([spt.hCat([uu,ux,ul]),\
                      spt.hCat([xu,xx,xl]),\
                      spt.hCat([lu,lx,ll])])
        
        self.projector = lin.factorized(M.tocsc())
        
    def contrastProjector(self,P,uTrue,xTrue):
        '''subroutine to actually do the projection '''
        n = self.fwd.N
        m = self.fwd.getXSize()

        # Aa = sparse.eye(n+m,n+m);
        # Bb = -sparse.eye(n+m,n+m);
        eAbs = 1e-5;
        eRel = 1e-5;

        u = np.zeros(n);
        ut = self.us;
        ud = np.zeros(n);

        x = np.zeros(m);
        xt = self.X;
        xd = np.zeros(m);

        z = np.zeros(n+m);

        TT = sparse.spdiags(self.s*self.fwd.p2x*P,0,m,m)*self.fwd.x2u.T
        
        uu = self.rho*TT.T.conj()*TT + self.xi*sparse.eye(n,n);
        ux = -self.rho*TT.T.conj()
        xu = -self.rho*TT
        xx = self.rho*sparse.eye(m,m) + self.xi*sparse.eye(m,m)
        
        R = spt.vCat([spt.hCat([uu, ux]), spt.hCat([xu, xx])])
        g = lin.factorized(R.tocsc())
        
        ePri = 0.0
        eDua = 0.0
        rErr = np.ones(n+m)
        sErr = np.ones(n+m)
        # gap = np.zeros(100)
        
        iter = 1
        while (iter<50): # & (np.linalg.norm(rErr)>ePri) & (np.linalg.norm(sErr)>eDua):
            ''' inner loop to solve the projection '''
            iter += 1
            rhs = np.concatenate((self.fwd.Ms.T*self.uHat + self.xi*(ut-ud),\
                                  self.xi*(xt-xd),\
                                  np.zeros(n)))
            updt = self.projector(rhs)
            u = updt[:n]
            x = updt[n:(n+m)]
            
            rhs = np.concatenate((self.xi*(u+ud) + self.rho*TT.T.conj()*(self.Z-TT*self.ub),\
                                  self.xi*(x+xd) - self.rho*(self.Z-TT*self.ub)))
            
            zold = z;
            z = g(rhs);
            ut = z[:n]
            xt = z[n:]        
    
            ud = ud + u-ut;
            xd = xd + x-xt;
            gap = np.linalg.norm(np.concatenate((u,x))-np.concatenate((ut,xt)))
            print 'Gap at iter ' + repr(iter) + ' ' + repr(gap)
            
            print 'uErr ' + repr(np.linalg.norm(u-uTrue)/np.linalg.norm(uTrue))
            print 'xErr ' + repr(np.linalg.norm(x-xTrue)/np.linalg.norm(xTrue))
            
            sErr = -self.rho*(z-zold);
            rErr = np.concatenate((u,x)) - z;
    
            ePri = np.sqrt(2*n)*eAbs + eRel*max(np.linalg.norm(np.concatenate((ut,xt))),\
                                                np.linalg.norm(-1.0*z))
                                                
            eDua = np.sqrt(2*n)*eAbs + eRel*np.linalg.norm(np.concatenate((ud,xd))*self.rho);
            # if (np.linalg.norm(rErr)<ePri) & (np.linalg.norm(sErr)<eDua):
            #    break
    
        print 'inner iters = ' + repr(iter)
#        plt.figure(1)
#        plt.plot(range(100),gap)
#        plt.show()
        return u,x
        

    def aggregatorSemiParallel(self,S, comm):
        ''' Do the aggregation step in parallel whoop! '''
        N = np.size(S)
        n = S[0].fwd.nRx*S[0].fwd.nRy
        
        uL = sparse.lil_matrix((n,n),dtype='complex128')
        bL = np.zeros(n,dtype='complex128')
        
#        U = np.zeros((n,N),dtype='complex128')
#        Q = np.zeros((n,N),dtype='complex128')
        nX = self.fwd.getXSize()
        for L in S:
            # for ix in range(N):
#            s = S[ix].s
            # print s
#            U[:,ix] = s*S[ix].fwd.Md*(S[ix].ub + S[ix].us)
#            Q[:,ix] = S[ix].X + S[ix].Z
            # print L.fwd.x2u.shape
            
            M = L.s*(sparse.spdiags(L.fwd.x2u.T*(L.us+L.ub),0,nX,nX))*self.fwd.p2x
            uL += M.T.conj()*M
            bL += M.T.conj()*(L.X + L.Z)
            
        U = sparse.lil_matrix((n,n),dtype='complex128')
        B = np.zeros(n,dtype='complex128')
        
        U = comm.allreduce(uL,U,op=MPI.SUM)
        B = comm.allreduce(bL,B,op=MPI.SUM)
        
        P = lin.spsolve(U,B)
        # print num
        # print den
        
        P = P.real
        
        # print P[1]
        P = np.maximum(P,0)
        # print self.upperBound
        P = np.minimum(P,self.upperBound)
        # print P[1]
               
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
        skt = self.uHat
        
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
