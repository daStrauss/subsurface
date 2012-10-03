'''
Created on Jun 12, 2012

@author: dstrauss
'''

from model import fwd
from scipy import sparse
import sparseTools as spt
import numpy as np
from scipy.sparse import linalg as lin
import scipy.io as spio

class solver(fwd):
    
    def makeGradOper(self):
        hz2ex = sparse.kron(sparse.eye(self.nx+1,self.nx+1), self.po*self.d2);
        hz2ey = sparse.kron(-self.po*self.d2, sparse.eye(self.nx+1,self.nx+1));

        ex2hz = sparse.kron(sparse.eye(self.nx+1,self.nx+1),-self.ph*self.d1);
        ey2hz = sparse.kron(self.ph*self.d1, sparse.eye(self.nx+1,self.nx+1));

        n = self.nx*(self.nx+1);
        N = (self.nx+1)**2;
        
        self.nabla2 = spt.vCat([spt.hCat([sparse.coo_matrix((n,n)), sparse.coo_matrix((n,n)), hz2ex]), \
                                spt.hCat([sparse.coo_matrix((n,n)), sparse.coo_matrix((n,n)), hz2ey]), \
                                spt.hCat([ex2hz, ey2hz, sparse.coo_matrix((N,N))]) ])
                              
    
                              
    def setmats(self, eHSr, sHS, div):
        """This quick routine starts the process to
        setup the location of the background halfspace 
        It is also true that this set the space for the ON GRID locations
        """
        self.eHS = eHSr*self.epso # half space permitivity
        self.sHS = sHS # half space conductivity
        self.div = div # hslf space dividing line
        self.kfree = 2*np.pi/self.l; # define k in free space
        self.kHS = np.sqrt(self.muo*self.eHS*(self.w**2) \
                         + 1j*self.w*self.muo*self.sHS); # k in subsurface
                         
        # ok. lets lose a little flexibility
        sigX = np.zeros((self.nx+1, self.nx))
        sigX[:,:(div+1)] = self.sHS
        
        sigY = np.zeros((self.nx,self.nx+1))
        sigY[:,:(div+2)] = self.sHS
        
        # and keep this as zeros- reason being, it is going to be much more effective, and lingustically simple
        # to not have to deal with doing things other than Md*u
        sigZ = np.zeros((self.nx+1,self.nx+1))
        
        # self.epsmap = [self.epso*np.ones((self.nx,self.ny)), self.epso*np.ones((self.nx,self.ny))]
        
        self.sigmap = np.concatenate((sigX.flatten(), sigY.flatten(), sigZ.flatten()))
        # and duplicate
        self.sigmap = [self.sigmap, self.sigmap.copy()]
        
        
        self.N = (self.ny+1)*self.nx + (self.nx+1)*self.nx + (self.nx+1)*(self.ny+1)
        self.sol = [np.zeros((self.N,1)), np.zeros((self.N,1))]

    
    def parseFields(self, u):
        ''' given a full space vector, u, it parses it into the three objects and returns them 
        [hz, ex, ey] -- could also be used for the sigmaps if wanted. a permutation, sure, but 
         '''
        hi = (self.nx+1)*(self.ny)

        ex = u[:hi]
        ex = ex.reshape(self.nx+1,self.ny)
        
        hj = hi + (self.nx)*(self.ny+1)
        ey = u[hi:hj]
        ey = ey.reshape(self.nx,self.ny+1)
        
        hz = u[hj:]
        hz = hz.reshape(self.nx+1,self.ny+1)
        
        return [hz,ex,ey]
        
        
            
            
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
        nLocal = (self.nx+1)*(self.ny+1)
        self.Ms = sparse.dok_matrix((nLocal,idx.size))
        
        # there's probably a more pythonic way here ---
        for i in range(sum(oprx)):
            self.Ms[idx[i],i] = 1.0
        self.Ms = self.Ms.tocsc()
        self.Ms = self.Ms.T
        
        z = self.Ms.shape
        #super dimension of the ex,ey
        n = (self.nx+1)*self.nx
        self.Ms = spt.hCat([sparse.coo_matrix((z[0],n)), sparse.coo_matrix((z[0],n)), self.Ms])
        
    def getk(self, idx):
        ''' method to return the diagonal for the self.nabla2 operator '''
        N = (self.nx+1)*self.ny + (self.ny+1)*self.nx
        n = (self.nx+1)*(self.ny+1)
        dia = np.append(1j*self.w*self.epso*np.ones(N), \
                        -1j*self.w*self.muo*np.ones(n)) \
                        - self.sigmap[idx] 
        return sparse.spdiags(dia, 0, self.N, self.N)
    
    def setCTRX(self):
        '''produce the operators p2x and x2u so that I can do the contrastX solves'''
        
        
        self.p2x = spt.vCat([sparse.eye(self.nRx*self.nRy,self.nRx*self.nRy), \
                             sparse.eye(self.nRx*self.nRy,self.nRx*self.nRy)])
        self.p2x = self.p2x.tocsc()
#        print self.p2x.shape
        
        N = (self.nx+1)*(self.ny+1)
        Nx = (self.nx+1)*self.ny
        Ny = self.nx*(self.ny+1)
        NR = self.nRx*self.nRy 
        
        lMd = self.Md.tolil()
        Mdx = lMd[:,:Nx].T.tocoo()
        Mdy = lMd[:,Nx:(Nx+Ny)].T.tocoo()
#        print 'sliced'
#        assert self.nRx*self.nRy == 
        self.x2u = spt.vCat([spt.hCat([Mdx, sparse.coo_matrix((Nx,NR))]),\
                             spt.hCat([sparse.coo_matrix((Ny,NR)), Mdy]),\
                             sparse.coo_matrix((N, 2*NR)) ])  #oh bugger, problems. 
#        print 'assembled'
        self.x2u = self.x2u.tocsc()
#        print 'converted'
#        self.x2u = self.x2u
#        print self.x2u.shape
        
    def getXSize(self):
        ''' return the proper size of X so that the optimizatoin routine can work its magic '''
        return 2*self.nRx*self.nRy
        
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
        Mdx = sparse.dok_matrix(((self.nx+1)*self.nx,idx.size))
        
        for i in range(idx.size):
            Mdx[idx[i],i]=1.0
            
        opry = np.zeros((self.nx,self.ny+1),dtype='bool')
        # might have to revise this later.
        opry[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1
        
        idx = np.arange((self.nx+1)*self.nx)
        opry = opry.flatten()
        idx = idx[opry]
        Mdy = sparse.dok_matrix(((self.nx+1)*self.nx,idx.size))
        
        for i in range(idx.size):
            Mdy[idx[i],i] = 1.0
        # only do one conversion to csc    
        N = (self.nx+1)*(self.nx+1)
        assert Mdx.shape == Mdy.shape
        
        
        
#        assert self.nRx*self.nRy == 
        self.Md = spt.vCat([Mdx, Mdy, sparse.coo_matrix((N,Mdx.shape[1]))]) #oh bugger, problems. 
        self.Md = self.Md.tocsc()
        self.Md = self.Md.T
#                              
# no need to overload the fwd_solve method -- its the same!
#    def fwd_solve(self, ind):
#        """ Does the clean solve, prints a figure """
#        # self.sol[ind] = self.rhs.copy();
#        # b = self.rhs.copy().flatten();
#        self.sol[ind] = lin.spsolve(sparse.csc_matrix(self.nabla2+ self.getk(ind)), \
#                                    self.rhs.flatten())
        
        # umfpack.linsolv((self.nabla2 + self.getk(ind)), self.sol[ind])
        # self.sol[ind] = np.array(self.sol[ind])

    def pointSource(self,x,y):
        ''' Make a magnetic field point source ''' 
        rhsz = np.zeros((self.nx+1,self.ny+1))
        rhsz[x,y] = -1
        
        rhsx = np.zeros((self.nx+1,self.ny))
        rhsy = np.zeros((self.nx,self.ny+1))

        self.rhs = np.concatenate((rhsx.flatten(), rhsy.flatten(), rhsz.flatten()))

                              
    def planeWave(self):
        ''' Make the plane wave illumination for the given frequency and inc angle '''

        #    The assumption is that the Ez and materials are co
        #    located. Since epsilon(50) => in the half space, epsilon(51) =>
        #    is not, the actual zero boundary must be between them, or on
        #    that y boundary.
        #Mapping is y,x because of how matlab treats these. annoying.
        #
        #%fully half grid
        x = np.arange(1,1+self.nx,dtype='float64')*self.dx
        cntr = (self.div+1+0.5)*self.dx
        Yhh,Xhh = np.meshgrid(np.append(0.0, x) + (self.dx/2) - cntr, \
                              np.append(0.0, x) + (self.dx/2)); 
                               
        # y - half grid
        Yyh, Xyh = np.meshgrid(np.append(0.0, x) + (self.dx/2) - cntr, \
                               x);
        
        #x - half grid
        Yxh,Xxh = np.meshgrid(x-cntr, \
                              np.append(0.0,x) + (self.dx/2));
  
        ni = 1.0
        nt = np.sqrt((self.eHS-self.sHS/(1j*self.w))*self.muo)/np.sqrt(self.epso*self.muo)
        #incident angle
        tht = np.arcsin(ni*np.sin(self.incAng)/nt)

        #Create the coefficients to specify the space. 
        kinc = -np.array([np.sin(self.incAng), np.cos(self.incAng)])
        ktx = -np.array([np.sin(tht), np.cos(tht)])
        kFS = 1j*np.sqrt(self.muo*self.epso*self.w**2)
        etaF = np.sqrt(self.muo/self.epso)

        kHS = 1j*np.sqrt(self.muo*self.eHS*self.w**2 + 1j*self.w*self.muo*self.sHS)
        etaH = np.sqrt(self.muo/(self.eHS+1j*self.sHS/self.w))
        
        rTM = (nt*np.cos(self.incAng) - ni*np.cos(tht))/(ni*np.cos(tht) + nt*np.cos(self.incAng));
        tTM = (2*ni*np.cos(self.incAng))/(ni*np.cos(tht) + nt*np.cos(self.incAng));
        
        # make space selector matrices
        ths = np.zeros((self.nx+1,self.ny+1), dtype='bool')
        # res.sigHz = zeros(nx+1,nx+1);
        ths[Yhh<=0] = 1
        #res.sigHz(Yhh<=0) = sHS;
        
        thsy = np.zeros((self.nx,self.nx+1), dtype='bool')
        #  res.sigEy = zeros(nx,nx+1);
        thsy[Yyh<=0] = 1
        #res.sigEy(Yyh<=0) = sHS;
        
        thsx = np.zeros((self.nx+1,self.nx), dtype='bool')
        #  res.sigEx = zeros(nx+1,nx);
        thsx[Yxh<=0] = 1
        #  res.sigEx(Yxh<=0) = sHS;
        
        #Make the background fields
        Hzinc = np.zeros((self.nx+1,self.ny+1), dtype='complex128');
        Hzinc[~ths] = (1/etaF)*np.exp(kFS*(Xhh[~ths]*kinc[0] + Yhh[~ths]*kinc[1])) + \
                  rTM*(1/etaF)*np.exp(kFS*(Xhh[~ths]*kinc[0] - Yhh[~ths]*kinc[1]))

        Hzinc[ths] = tTM*(1/etaH)*np.exp(kHS*(Xhh[ths]*ktx[0] + Yhh[ths]*ktx[1]))
        
        Exinc = np.zeros((self.nx+1,self.ny),dtype='complex128');
        Exinc[~thsx] = (kinc[1]*np.exp(kFS*(Xxh[~thsx]*kinc[0] + Yxh[~thsx]*kinc[1]))) + \
                  rTM*(-kinc[1]*np.exp(kFS*(Xxh[~thsx]*kinc[0] - Yxh[~thsx]*kinc[1])))

        Exinc[thsx] = tTM*(ktx[1]*np.exp(kHS*(Xxh[thsx]*ktx[0] + Yxh[thsx]*ktx[1])))
        Exinc = -   Exinc
        
        Eyinc = np.zeros((self.nx,self.ny+1), dtype='complex128')
        Eyinc[~thsy] = (kinc[0]*np.exp(kFS*(Xyh[~thsy]*kinc[0] + Yyh[~thsy]*kinc[1]))) + \
                   rTM*(kinc[0]*np.exp(kFS*(Xyh[~thsy]*kinc[0] - Yyh[~thsy]*kinc[1])))
                   
        Eyinc[thsy] = tTM*(ktx[0]*np.exp(kHS*(Xyh[thsy]*ktx[0] + Yyh[thsy]*ktx[1])))
        
        
        bnd = [Exinc, Eyinc,Hzinc]
#        Eyinc = -Eyinc;
        
#        fldz.Hz = Hzinc;
#        fldz.Ex = Exinc;
#        fldz.Ey = Eyinc;
  
#        fldz.sEx = reshape(hz2ex*Hzinc(:), nx+1,nx)./(-1i*w*epso+res.sigEx);
#        fldz.sEy = reshape(hz2ey*Hzinc(:), nx,nx+1)./(-1i*w*epso+res.sigEy);
#        fldz.sHz = reshape(ex2hz*Exinc(:) + ey2hz*Eyinc(:), nx+1,nx+1)/(1i*w*muo);

        Msrcz = np.zeros((self.nx+1,self.ny+1),dtype='complex128')
        Jsrcx = np.zeros((self.nx+1,self.ny), dtype='complex128')
        Jsrcy = np.zeros((self.nx,self.ny+1), dtype='complex128')
        
        instep = 3+self.npml;
        xl = instep-1
        xr = self.nx-instep-1
        yb = instep-1
        yt = self.ny-instep-1

#        5.48a
        Jsrcx[(xl+1):(xr+1),yb] += (Hzinc[(xl+1):(xr+1),yb]/self.dx)
#        5.49a
        Jsrcx[(xl+1):(xr+1),yt] += ((-1.0)*Hzinc[(xl+1):(xr+1),yt+1]/self.dx)

#        5.52a
        Jsrcy[xl,(yb+1):(yt+1)] += ((-1.0)*Hzinc[xl,(yb+1):(yt+1)]/self.dx)
#        5.53a
        Jsrcy[xr,(yb+1):(yt+1)] += (Hzinc[xr+1,(yb+1):(yt+1)]/self.dx)

#        5.54a
        Msrcz[(xl+1):(xr+1),yb] += ((-1.0)*Exinc[(xl+1):(xr+1),yb]/self.dx)
        
#        5.55a
        Msrcz[(xl+1):(xr+1),yt+1] += (Exinc[(xl+1):(xr+1),yt]/self.dx)
  
#        5.58a
        Msrcz[xl, (yb+1):(yt+1)] += (Eyinc[xl, (yb+1):(yt+1)]/self.dx)

#        5.59a
        Msrcz[xr+1,(yb+1):(yt+1)] += ((-1.0)*Eyinc[xr, (yb+1):(yt+1)]/self.dx)
        
        
        self.rhs = np.concatenate((Jsrcx.flatten(), Jsrcy.flatten(), Msrcz.flatten()))
        return bnd
    
    def getS(self):
        '''return the parameter for scaling within the Md*p part. '''
        return -1.0+0j
