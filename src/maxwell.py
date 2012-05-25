'''
Created on May 23, 2012

@author: dstrauss
'''
import numpy as np
# import scipy as sp
from scipy import sparse
from scipy.sparse import linalg as lin
import matplotlib.pyplot as plt
import scipy.special as spec
import pickle

class pmlList(object):
    freq = np.ndarray(0)
    pmc = np.ndarray(0)
            

class twoDim(object):
    ''' 
    A two dimensional Maxwell equation simulator class for finite difference frequency domain work
    '''
    
    # define and keep universal constants
    epso = 8.854e-12
    muo = 4.0*np.pi*1e-7
    c = np.sqrt(1/(muo*epso))
    
    # create a few holders for some important things
    sol = ['','','']
    # A = ['','','']
    Q = ['','','']
    
    def __init__(self, freq):
        self.w = 2*np.pi*freq
        self.f = freq
        self.l = self.c/self.f
        print 'building a f = %d object' %freq
        
    def setspace(self, nx,ny,dx,dy):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.N = nx*ny # super size of space
        self.dx = dx # delta x
        self.dy = dy # delta y
        # interestingly, I want to not include zeros((q),1) because that gives
        # the wrong size arrays
        self.rhs = np.zeros(self.N, dtype='complex128')
        # np.zeros((self.nx,self.nx),complex);
        self.npml = min(10,round((nx+2.0)/10)) # size of pml borders
        
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

        self.epsmap = [self.epso*np.ones((self.nx,self.ny)), self.epso*np.ones((self.nx,self.ny))]
        self.sigmap = [np.zeros((self.nx,self.ny)), np.zeros((self.nx,self.ny))]

        for x in range(2):
            self.epsmap[x][:,0:div] = self.eHS
            self.sigmap[x][:,0:div] = self.sHS
        
    def getk(self, ind):
        """ This routine assembles a diagonal matrix with the materials indexed by ind
        """
        kl = (self.muo*self.epsmap[ind]*(self.w**2) + 1j*self.w*self.muo*self.sigmap[ind]   )
        
        return sparse.spdiags(kl.flatten(), [0], self.N, self.N)
        
    def prmz(self):
        """This routine might be nice to have it spit out a basic
        parameters in the structure """
        print('Size ({0}, {1}); gridsize ({2}, {3})'.format(self.nx,self.ny,self.dx,self.dy))
    
    def getPMLparm(self):
        ''' Returns the proper PML parameter either from a library if it has already been calculated,
        or it is able to explicitly calculate it otherwise
        '''
        try:
            lkt = pickle.load(open('pmlLib.p', 'rb'))
        except:
            lkt = pmlList()
            
        if np.any(lkt.freq==self.f):
            return lkt.pmc[lkt.freq==self.f]
        else:
            pmc = findBestAng(self.f)
            lkt.freq = np.append(lkt.freq, self.f)
            lkt.pmc = np.append(lkt.pmc, pmc)
            pickle.dump(lkt, open('pmlLib.p', 'wb'))
            return pmc
            
    def setOperators(self):
        ''' An wraper for makeFinteDiferences that will select the right PML parameters
        '''
        self.makeFD(self.getPMLparm())
            
            
        
    def makeFD(self, pmc):
        """A routine that will define the following operators:
        A, ox, oy"""
        # h = (self.dx)*np.mat(np.ones((self.nx+1,1)));
        hi = np.zeros(self.nx+1);
    
        hi[:self.npml] = (pmc)*np.linspace(1,0,self.npml);
        hi[(self.nx-self.npml+1):(self.nx+1)] = \
                  (pmc)*np.linspace(0,1,self.npml);

        h = np.vectorize(complex)(self.dx*np.ones(self.nx+1), hi);


        opr = np.ones((2,self.nx+1))
        opr[0,:] = -1
        
        self.d1 = sparse.spdiags(opr, [-1, 0], self.nx+1, self.nx);
        self.d2 = sparse.spdiags(opr, [0,  1], self.nx, self.nx+1);

        # averaging operator for getting to the in-between grid points
        avg_n_c = sparse.spdiags(0.5*np.ones((2,self.nx+1)), [0, 1], self.nx,self.nx+1);

        #creating the half and non-half space streatchers.
        ph = sparse.spdiags(1/h, 0, self.nx+1,self.nx+1);
        po = sparse.spdiags(1/(avg_n_c*h), 0, self.nx,self.nx);

        oper1d = po*self.d2*ph*self.d1;
        # ok. building nabla2 for the TE only
        self.nabla2 = sparse.kron(oper1d,sparse.eye(self.nx,self.nx)) \
                 + sparse.kron(sparse.eye(self.nx,self.nx), oper1d);

        
    def point_source(self, x,y):
        """ A routine to add a point source at the grid loc (x,y) """
        self.rhs.reshape(self.nx,self.ny)[x,y] = -1.0

    def te_pw(self, thi):
        """ A routine to add a te planewave at angle as spec'd """
        instep = 3+(self.nx+1)/10;
        # mdpt = nx/2; # should replace by div
        x = np.arange(1,self.nx+1)*self.dx
        # The assumption is that the Ez and materials are co
        # located. Since epsilon(50) => in the half space, epsilon(51) =>
        # is not, the actual zero boundary must be between them, or on
        # that y boundary.
        # Mapping is y,x because of how matlab treats these. annoying.
        # (
        Y,X = np.meshgrid(x-(self.div+0.5)*self.dx, x); # I do things backward
        Yyh, Xyh = np.meshgrid(np.append(0.0,x)+(self.dx/2)\
                               - (self.div+0.5)*self.dx,x);
        Yxh, Xxh = np.meshgrid(x-(self.div+0.5)*self.dx, \
                               np.append(0.0,x)+(self.dx/2));
  
        ni = 1;
        nt = np.sqrt((self.eHS-self.sHS/(1j*self.w))\
                     *self.muo)/np.sqrt(self.epso*self.muo);

        # transmitted angle angle
        # thi = 45*np.pi/180 # taken as input argument.
        tht = np.arcsin(ni*np.sin(thi)/nt);

        # create the coefficients to specify the space. 
        kinc = -1*np.array([np.sin(thi), np.cos(thi)])
        ktx = -1*np.array([np.sin(tht), np.cos(tht)])
        kFS = 1j*self.kfree;
        kHS = 1j*self.kHS;
        etaF = np.sqrt(self.muo/self.epso);
        #   % [kFS kfs]
        etaH = np.sqrt(self.muo/(self.eHS+1j*self.sHS/self.w));
        rTE = (ni*np.cos(thi) - nt*np.cos(tht))/(ni*np.cos(thi) + nt*np.cos(tht));
        tTE = (2*ni*np.cos(thi))/(ni*np.cos(thi) + nt*np.cos(tht));

        #   % [yy xx] = meshgrid(zeY, zeX);
        #   % Make a selector for the half space.
        ths = np.zeros([self.nx,self.ny])
        ths[Y<0] = 1
        ths = ths.astype(bool)
        #   size(ths);
  
        thsy = np.zeros([self.nx,self.nx+1])
        thsy[Yyh<0] = 1
        thsy = thsy.astype(bool)

        thsx = np.zeros([self.nx+1,self.nx])
        thsx[Yxh<0] = 1
        thsx = thsx.astype(bool)
  
        # % kHS
        # %  size(thsy)
        Ezinc = np.zeros([self.nx,self.ny], dtype='complex128')
        Ezinc[~ths] = np.exp(kFS*(X[~ths]*kinc[0] + Y[~ths]*kinc[1])) + \
              rTE*np.exp(kFS*(X[~ths]*kinc[0] - Y[~ths]*kinc[1]))
  
        Ezinc[ths] = tTE*np.exp(kHS*(X[ths]*ktx[0] + Y[ths]*ktx[1]))
  
        Hxinc = np.zeros([self.nx,self.nx+1],dtype='complex128');
        Hxinc[~thsy] = (1/etaF)*(kinc[1]*np.exp(kFS*(Xyh[~thsy]*kinc[0] + Yyh[~thsy]*kinc[1]))) \
               +rTE*(1/etaF)*(-kinc[1]*np.exp(kFS*(Xyh[~thsy]*kinc[0] - Yyh[~thsy]*kinc[1])))
  
        Hxinc[thsy] = tTE*(1/etaH)*(ktx[1]*np.exp(kHS*(Xyh[thsy]*ktx[0] + Yyh[thsy]*ktx[1])));
        #   Hxinc = Hxinc;
  
        Hyinc = np.zeros([self.nx+1,self.nx], dtype='complex128');
        Hyinc[~thsx] = (1/etaF)*(kinc[0]*np.exp(kFS*(Xxh[~thsx]*kinc[0] + \
                                             Yxh[~thsx]*kinc[1]))) + \
               rTE*(1/etaF)*(kinc[0]*np.exp(kFS*(Xxh[~thsx]*kinc[0] - \
                                                 Yxh[~thsx]*kinc[1])))

        Hyinc[thsx] = tTE*(1/etaH)*(ktx[0]*np.exp(kHS*(Xxh[thsx]*ktx[0] + \
                                            Yxh[thsx]*ktx[1])))
        Hyinc = -Hyinc;
        
#        plt.figure(11)
#        plt.subplot(1,3,1)
#        plt.imshow(Ezinc.real)
#        plt.title('Real Ez inc')
#        plt.colorbar()
#        
#        plt.subplot(1,3,2)
#        plt.imshow(Hxinc.real)
#        plt.title('Real Hx inc')
#        plt.colorbar()
#
#        plt.subplot(1,3,3)
#        plt.imshow(Hyinc.imag)
#        plt.title('Imag Hy inc')
#        plt.colorbar()

        xl = instep-1; xr = self.nx-instep-1;
        yb = instep-1; yt = self.ny-instep-1;
        Jsrcz = np.zeros([self.nx,self.ny],dtype='complex128');
        Msrcx = np.zeros([self.nx,self.ny+1], dtype='complex128');
        Msrcy = np.zeros([self.nx+1,self.ny], dtype='complex128');

        Jsrcz[xl,yb:yt] =    Jsrcz[xl,yb:yt] + (1)*(Hyinc[xl,yb:yt]/self.dx);
        Jsrcz[xr,yb:yt] =    Jsrcz[xr,yb:yt] - (1)*(Hyinc[xr+1,yb:yt]/self.dx);
        Jsrcz[xl:xr,yb] =    Jsrcz[xl:xr,yb] - (1)*(Hxinc[xl:xr,yb]/self.dy);
        Jsrcz[xl:xr,yt] =    Jsrcz[xl:xr,yt] + (1)*(Hxinc[xl:xr,yt+1]/self.dy);

        Msrcx[xl:xr,yb]   =  (1)*(Ezinc[xl:xr,yb]/self.dy);
        Msrcx[xl:xr,yt+1] = -(1)*(Ezinc[xl:xr,yt]/self.dy);
    
        Msrcy[xl,   yb:yt] = -(1)*(Ezinc[xl,yb:yt]/self.dx);
        Msrcy[xr+1, yb:yt] =  (1)*(Ezinc[xr,yb:yt]/self.dx);
        
#        aliJ = self.nx*self.ny
#        aliMx = self.nx*(self.ny+1)
#        aliMy = (self.nx+1)*self.ny
        print Msrcx.shape
        print Msrcy.shape
        
        # note -- this is different from the MATLAB version: due to the row-ordering of arrays
        pw = (1j*self.w*self.muo)*Jsrcz.flatten() - \
              sparse.kron(sparse.eye(self.nx,self.nx), self.d2/self.dx)*Msrcx.flatten() + \
              sparse.kron(self.d2/self.dx, sparse.eye(self.ny,self.ny))*Msrcy.flatten()

        self.rhs = self.rhs + pw.flatten();

    def fwd_solve(self, ind):
        """ Does the clean solve, prints a figure """
        # self.sol[ind] = self.rhs.copy();
        # b = self.rhs.copy().flatten();
        self.Q[ind] = lin.factorized(sparse.csc_matrix(self.nabla2 + self.getk(ind)))
        
        self.sol[ind] = self.Q[ind](self.rhs.flatten())
        # umfpack.linsolv((self.nabla2 + self.getk(ind)), self.sol[ind])
        # self.sol[ind] = np.array(self.sol[ind])
        self.sol[ind] = self.sol[ind].reshape(self.nx,self.ny)
        
    def setMs(self, nSensors=10):
        '''Tell me the number of sensors, and I will distribute them equally across the surface
        '''
        indx = np.round(np.linspace(self.npml+10,self.nx-self.npml-10, nSensors)).astype(int);
        oprx = np.zeros((self.nx,self.ny),dtype='bool')
        
        oprx[self.div+1,indx] = 1;
        
        idx = np.arange(self.N)
        oprx = oprx.flatten()
        idx = idx[oprx]
        
        self.Ms = sparse.dok_matrix((self.N,idx.size), dtype='bool')
        
        for i in range(sum(oprx)):
            self.Ms[idx[i],i] = 1
        self.Ms = self.Ms.tocsc()
        self.Ms = self.Ms.T
        
        
    def setMd(self, xrng, yrng):
        '''Tell me the xrange and the yrange and I'll make selector
        '''
        oprx = np.zeros((self.nx,self.ny),dtype='bool')
        oprx[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1
        self.nRx = xrng[1]-xrng[0] + 1
        self.nRy = yrng[1]-yrng[0] + 1
        
        idx = np.arange(self.N)
        oprx = oprx.flatten()
        idx = idx[oprx]
        self.Md = sparse.dok_matrix((self.N,idx.size), dtype = 'bool')
        
        for i in range(idx.size):
            self.Md[idx[i],i]=1
        self.Md = self.Md.tocsc()
        self.Md = self.Md.T
        

            
    def plotSol(self,ind):
        pass
    
def findBestAng(freq):
    '''for a particular frequency, find the best PML complex angle 
    '''
    # set some internal parameters
    nx = 149
    ny = nx
    dx = 5.0
    dy = dx
    
    # this means free space
    eHS = 1.0
    sHS = 0
    
    epso = 8.854e-12
    muo = 4.0*np.pi*1e-7
    
    x = np.array(range(1,nx+1))*dx - dx*(nx/2)
    Y,X = np.meshgrid(x, x);
    dist = np.sqrt(Y**2+X**2) 
    
    
    
    c = 1/np.sqrt(muo*epso)
    
    k = (2*np.pi)/(c/freq)
    bbx = dx*dx*(1j/4)*spec.hankel1(0,k*dist)
    mdpt = nx/2-1
    # the center contains a singularity, but we'll  be forgiving and finite
    bbx[mdpt,mdpt] = 0.25*bbx[mdpt+1,mdpt] + 0.25*bbx[mdpt-1,mdpt] + 0.25*bbx[mdpt,mdpt+1] + 0.25*bbx[mdpt,mdpt-1]
    
    sze = bbx.shape
    mask = np.ones(sze)
    mask[:15,:] =0
    mask[(nx-15):(nx-1),:] = 0
    mask[:,:15]=0
    mask[:,(nx-15):(nx-1)] =0
    
    lo = 0
    hi = 5
    
    # do two loops over different ranges to get a more precise estimate
    for tune in range(2):
        localError = np.zeros(50)
        angChoice = np.logspace(lo,hi,50)
    
        for i in range(50):    
            #     create a single object
            bce = twoDim(freq)
            bce.setspace(nx,ny,dx,dy)
            bce.setmats(eHS, sHS, nx/2)
    
            # pmc = 22.229964825261945
            bce.makeFD(angChoice[i])
            bce.point_source(nx/2 -1, ny/2 - 1)
            bce.fwd_solve(0)
    
            localError[i] = np.linalg.norm((bce.sol[0] - bbx)*mask,'fro')
            print 'ang = ' + repr(angChoice[i]) + ' local error ' + repr(localError[i])
    
        minIdx = np.argmin(localError)
        lo = np.log10(angChoice[max(0,minIdx-1)])
        hi = np.log10(angChoice[min(50,minIdx+1)])
        
    bce = twoDim(freq)
    bce.setspace(nx,ny,dx,dy)
    bce.setmats(eHS, sHS, nx/2)
    
            # pmc = 22.229964825261945
    bce.makeFD(angChoice[minIdx])
    bce.point_source(nx/2 -1, ny/2 - 1)
    bce.fwd_solve(0)
    
    #do some representative plots
#    plt.figure(1)
#    plt.plot(angChoice, localError)
#    # do some plotting
#    plt.figure(13)
#    plt.subplot(1,2,1)
#    plt.imshow((bce.sol[0].real*mask))
#    plt.colorbar()
#    
#    plt.subplot(1,2,2)
#    plt.imshow((bbx.real*mask))
#    plt.colorbar()
#    
#    plt.figure(4)
#    plt.subplot(121)
#    plt.imshow((bce.sol[0]-bbx).real)
#    plt.colorbar()
#    plt.show()
    
    return angChoice[minIdx]
    
    