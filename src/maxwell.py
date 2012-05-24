'''
Created on May 23, 2012

@author: dstrauss
'''
import numpy as np
# import scipy as sp
from scipy import sparse
from scipy.sparse import linalg as lin

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
    A = ['','','']
    
    def __init__(self, freq):
        self.w = 2*np.pi*freq
        self.f = freq
        self.l = self.c/self.f
        print 'building a f = %d object' %freq
        
    def setspace(self, nx,ny,dx,dy):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.dx = dx # delta x
        self.dy = dy # delta y
        # interestingly, I want to not include zeros((q),1) because that gives
        # the wrong size arrays
        self.rhs = np.zeros((self.nx*self.ny), dtype='complex128')
        # np.zeros((self.nx,self.nx),complex);
        self.npml = min(10,round((nx+2)/10)) # size of pml borders
        
    def setmats(self, eHSr, sHS, div):
        """This quick routine starts the process to
        setup the location of the background halfspace 
        It is also true that this set the space for the ON GRID locations
        """
        self.eHS = eHSr*self.epso
        self.sHS = sHS
        self.div = div
        self.kfree = 2*np.pi/self.l;
        self.kHS = np.sqrt(self.muo*self.eHS*(self.w**2) \
                         + 1j*self.w*self.muo*self.sHS);

        self.epsmap = [self.epso*np.ones((self.nx,self.ny)), self.epso*np.ones((self.nx,self.ny))]
        self.sigmap = [np.zeros((self.nx,self.ny)), np.zeros((self.nx,self.ny))]

        for x in range(2):
            self.epsmap[x][:,0:div] = self.eHS
            self.sigmap[x][:,0:div] = self.sHS
        
    def getk(self, ind):
        """ This routine assembles a diagonal matrix with the materials indexed by ind
        """
        kl = (self.muo*self.epsmap[ind].T*(self.w**2) \
                     + 1j*self.w*self.muo*self.sigmap[ind].T)
        return sparse.spdiags((kl.reshape(self.nx*self.ny,1)).T, [0], self.nx*self.ny, self.nx*self.ny)
        
    def prmz(self):
        """This routine might be nice to have it spit out a basic
        parameters in the structure """
        print('Size ({0}, {1}); gridsize ({2}, {3})'.format(self.nx,self.ny,self.dx,self.dy))
        
    def make_operators(self, pmc):
        """A routine that will define the following operators:
        A, ox, oy"""
        # h = (self.dx)*np.mat(np.ones((self.nx+1,1)));
        hi = np.mat(np.zeros((self.nx+1,1)));
    
        hi[:self.npml] = (pmc)*np.linspace(1,0,self.npml).reshape(self.npml,1);
        hi[(self.nx-self.npml+1):(self.nx+1)] = \
                  (pmc)*np.linspace(0,1,self.npml).reshape(self.npml,1);

        h = np.vectorize(complex)(self.dx*np.ones((self.nx+1,1)), hi);

        opr = np.tile([-1,1], ((self.nx)+1,1));
        d1 = sparse.spdiags(opr.T, [-1, 0], self.nx+1, self.nx);
        d2 = sparse.spdiags(opr.T, [0,  1], self.nx, self.nx+1);
        self.d1 = d1 # hold on to these for later - they are 1d zero dirichlet operators
        self.d2 = d2 # hold on to these for later - 1d zdc oper
        
        avg_n_c = sparse.spdiags(np.tile([0.5, 0.5], (self.nx+1,1)).T, [0, 1], self.nx,self.nx+1);

        ph = sparse.spdiags(1/h.T, 0, self.nx+1,self.nx+1);
        po = sparse.spdiags(1/(avg_n_c*h).T, 0, self.nx,self.nx);

        oper1d = po*d2*ph*d1;
        oper2d = sparse.kron(oper1d,sparse.eye(self.nx,self.nx)) \
                 + sparse.kron(sparse.eye(self.nx,self.nx), oper1d);
        oper2d = oper2d.tocoo();

        # go ahead and keep it internally rather than return it.
        self.nabla2 = oper2d.copy() #spmatrix(oper2d.data,oper2d.row,oper2d.col)
        
    def point_source(self, x,y):
        """ A routine to add a point source at the grid loc (x,y) """
        self.rhs.reshape(self.nx,self.ny)[x,y] = -1.0

    def te_pw(self, thi):
        """ A routine to add a te planewave at angle as spec'd """
        instep = 3+(self.nx+1)/10;
        # mdpt = nx/2; # should replace by div
        x = np.array(range(1,self.nx+1))*self.dx
        # The assumption is that the Ez and materials are co
        # located. Since epsilon(50) => in the half space, epsilon(51) =>
        # is not, the actual zero boundary must be between them, or on
        # that y boundary.
        # Mapping is y,x because of how matlab treats these. annoying.
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
        
        # plt.figure(11)
        # plt.subplot(1,3,1)
        # plt.imshow(Ezinc.real)
        # plt.title('Real Ez inc')
        # plt.colorbar()
        
        # plt.subplot(1,3,2)
        # plt.imshow(Hxinc.real)
        # plt.title('Real Hx inc')
        # plt.colorbar()

        # plt.subplot(1,3,3)
        # plt.imshow(Hyinc.imag)
        # plt.title('Imag Hy inc')
        # plt.colorbar()

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
        
        aliJ = self.nx*self.ny
        aliMx = self.nx*(self.ny+1)
        aliMy = (self.nx+1)*self.ny
        pw = (1j*self.w*self.muo)*Jsrcz.T.reshape(aliJ,1) - \
              sparse.kron(self.d2/self.dx, \
                          sparse.eye(self.nx,self.nx))*Msrcx.T.reshape(aliMx,1) + \
              sparse.kron(sparse.eye(self.ny,self.ny), \
                          self.d2/self.dx)*Msrcy.T.reshape(aliMy,1)

        self.rhs = self.rhs + pw;

    def fwd_solve(self, ind):
        """ Does the clean solve, prints a figure """
        self.sol[ind] = self.rhs.copy();
        b = self.rhs.copy();
        self.Q[ind] = lin.factorized(self.nabla2 + self.getk(ind))
        
        self.sol[ind] = self.Q(b)
        # umfpack.linsolv((self.nabla2 + self.getk(ind)), self.sol[ind])
        self.sol[ind] = np.array(self.sol[ind])
        self.sol[ind] = self.sol[ind].reshape(self.nx,self.ny)
        
    def setSensorOps(self, nSensors):
        pass