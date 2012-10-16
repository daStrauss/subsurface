'''
Created on Oct 10, 2012

@author: dstrauss
'''

'''
Created on Jun 12, 2012

@author: dstrauss
'''

from model import fwd
import numpy as np
from scipy import sparse
import sparseTools as spt

def speye(n):
    return sparse.eye(n,n)

class solver(fwd):
    ''' class to implement the transverse electric mode, rather - the case where we have Ez only ''' 
    def setspace(self, nx,ny,nz,dx, dy, dz):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.nz = nz # number z points
        self.N = nx*ny*nz # super size of space
        self.dx = dx # delta x
        self.dy = dy # delta y
        self.dz = dz # delta z

        self.npml = min(10,round((nx+2.0)/10))
    
    def makeGradOper(self):
        ''' routine to make a big matrix for TE problems ex,ey,ez all incorporated,
        based on the petsc_cpx routine.''' 
        # a quick hint: pd2 == pdo
        pd1 = self.ph*self.d1
        pd2 = self.po*self.d2
        
        AA = sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd2))*\
        sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd1)) + \
        sparse.kron(pd2,speye((self.nx+1)*self.ny))*sparse.kron(pd1,speye((self.nx+1)*self.ny))

        AB = sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd2))*\
            sparse.kron(speye(self.nz),sparse.kron(pd1,speye(self.ny+1)))
            
        AC = sparse.kron(pd2,speye(self.ny*(self.nx+1)))*\
            sparse.kron(speye(self.nz+1),sparse.kron(pd1,speye(self.ny)))

        BA = sparse.kron(speye(self.nz),sparse.kron(pd2,speye(self.nx+1)))*\
            sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd1))
# chng
        BB = sparse.kron(speye(self.nz),sparse.kron(pd2,speye(self.nx+1)))*\
            sparse.kron(speye(self.nz),sparse.kron(pd1,speye(self.nx+1))) + \
            sparse.kron(pd2,speye((self.ny+1)*self.nx))*sparse.kron(pd1,speye((self.ny+1)*self.nx))
# chng
        BC = sparse.kron(pd2,speye((self.ny+1)*self.nx))*\
            sparse.kron(speye(self.nz+1),sparse.kron(speye(self.nx),pd1))
# chng
        CA = sparse.kron(speye(self.nz+1),sparse.kron(pd2,speye(self.ny)))*\
            sparse.kron(pd1,speye((self.nx+1)*self.ny))
# chng
        CB = sparse.kron(speye(self.nz+1),sparse.kron(speye(self.nx),pd2))*\
            sparse.kron(pd1,speye((self.ny+1)*self.nx))
# chng
        CC = sparse.kron(speye(self.nz+1),sparse.kron(speye(self.nx),pd2))*\
            sparse.kron(speye(self.nz+1),sparse.kron(speye(self.nx),pd1)) + \
            sparse.kron(speye(self.nz+1),sparse.kron(pd2,speye(self.ny)))*\
            sparse.kron(speye(self.nz+1),sparse.kron(pd1,speye(self.ny)))
            # chng
            
           # legacy - matlab ordering 
#    AA = sparse.kron(speye(self.nz),sparse.kron(pd2,speye(self.nx+1)))*\
#        sparse.kron(speye(self.nz),sparse.kron(pd1,speye(self.nx+1))) + \
#        sparse.kron(pd2,speye((self.nx+1)*self.ny))*sparse.kron(pd1,speye((self.nx+1)*self.ny))
#
#        AB = sparse.kron(speye(self.nz),sparse.kron(pd2,speye(self.nx+1)))*\
#            sparse.kron(speye(self.nz),sparse.kron(speye(self.ny+1), pd1))
#            
#        AC = sparse.kron(pd2,speye(self.ny*(self.nx+1)))*\
#            sparse.kron(speye(self.nz+1),sparse.kron(speye(self.ny), pd1))
#
#        BA = sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd2))*\
#            sparse.kron(speye(self.nz),sparse.kron(pd1,speye(self.nx+1)))
#
#        BB = sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd2))*\
#            sparse.kron(speye(self.nz),sparse.kron(speye(self.nx+1),pd1)) + \
#            sparse.kron(pd2,speye((self.ny+1)*self.nx))*sparse.kron(pd1,speye((self.ny+1)*self.ny))
#
#        BC = sparse.kron(pd2,speye((self.ny+1)*self.nx))*\
#            sparse.kron(speye(self.nz+1),sparse.kron(pd1,speye(self.nx)))
#
#        CA = sparse.kron(speye(self.nz+1),sparse.kron(speye(self.ny),pd2))*\
#            sparse.kron(pd1,speye((self.nx+1)*self.ny))
#
#        CB = sparse.kron(speye(self.nz+1),sparse.kron(pd2,speye(self.nx)))*\
#            sparse.kron(pd1,speye((self.ny+1)*self.nx))
#
#        CC = sparse.kron(speye(self.nz+1),sparse.kron(pd2,speye(self.nx)))*\
#            sparse.kron(speye(self.nz+1),sparse.kron(pd1,speye(self.nx))) + \
#            sparse.kron(speye(self.nz+1),sparse.kron(speye(self.ny),pd2))*\
#            sparse.kron(speye(self.nz+1),sparse.kron(speye(self.ny),pd1))

     
        self.nabla2  = spt.vCat([spt.hCat([AA, -AB, -AC]),\
                            spt.hCat([-BA, BB, -BC]), \
                            spt.hCat([-CA, -CB, CC])])

        
    
    def setmats(self, eHSr, sHS, div):
        """This routine does setup for the diagonal entries of the nabla matrix.
        Also included are interpolation operators to mapf rom half spaces to non-half spaces, maybe
            
        """
        self.eHS = eHSr*self.epso # half space permitivity
        self.sHS = sHS # half space conductivity
        self.div = div # hslf space dividing line
        self.kfree = 2*np.pi/self.l; # define k in free space
        self.kHS = np.sqrt(self.muo*self.eHS*(self.w**2) \
                         + 1j*self.w*self.muo*self.sHS); # k in subsurface
                         
        self.N = (self.nx+1)*self.ny*self.nz + \
                    self.nx*(self.ny+1)*self.nz + \
                    self.nx*self.ny*(self.nz+1)
        
        
        # epsilon is easy it is always going to be uniform
        self.epsmap = [self.epso*np.ones(self.N),\
                       self.epso*np.ones(self.N)]
        
        # sigma is not so straight forward
        
        sigX = np.zeros((self.nx+1,self.ny,self.nz))
        sigX[:,:(div+1),:] = self.sHS
        
        sigY = np.zeros((self.nx,self.ny+1,self.nz))
        sigY[:,:(div+2),:] = self.sHS
        
        sigZ = np.zeros((self.nx,self.ny,self.nz+1))
        sigZ[:,:(div+1),:] = self.sHS
        
        
        self.sigmap = np.concatenate((sigX.flatten(), sigY.flatten(), sigZ.flatten()))
        # and duplicate
        self.sigmap = [self.sigmap, self.sigmap.copy()]
        
        # this is the total number of unknown field values in the entire simulation space
        self.sol = [np.zeros((self.N,1)), np.zeros((self.N,1))]
        
            
    def getk(self, ind):
        """ This routine assembles a diagonal matrix with the materials indexed by ind
        """
        kl = (self.muo*self.epsmap[ind]*(self.w**2) + 1j*self.w*self.muo*self.sigmap[ind]   )
        
        return sparse.spdiags(kl.flatten(), 0, self.N, self.N)
    
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
        
        self.Ms = sparse.lil_matrix((self.N,idx.size))
        
        for i in range(sum(oprx)):
            self.Ms[idx[i],i] = 1.0
        self.Ms = self.Ms.tocsc()
        self.Ms = self.Ms.T
        
    def setCTRX(self):
        ''' create some operators to map back and forth between the x space and the u '''
#        self.p2x = sparse.eye(self.nRx*self.nRy,self.nRx*self.nRy)
#        # print self.p2x.shape
#        self.x2u = self.Md.T
#        # print self.x2u.shape
    
    def getXSize(self):
        ''' return the proper size of X so that the optimization routine can work its magic '''
        return self.nRx*self.nRy*self.nRz
        
    def setMd(self, xrng, yrng, zrng):
        '''Tell me the xrange and the yrange and I'll make selector'''
#        oprx = np.zeros((self.nx,self.ny),dtype='bool')
#        oprx[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1
#        self.nRx = xrng[1]-xrng[0]
#        self.nRy = yrng[1]-yrng[0]
#        
#        idx = np.arange(self.N)
#        oprx = oprx.flatten()
#        idx = idx[oprx]
#        self.Md = sparse.dok_matrix((self.N,idx.size))
#        
#        for i in range(idx.size):
#            self.Md[idx[i],i]=1.0
#        self.Md = self.Md.tocsc()
#        self.Md = self.Md.T
        
    def parseFields(self,u):
        ''' Method to return the field in its square form'''
        hi = (self.nx+1)*(self.ny)*self.nz

        ex = u[:hi]
        ex = ex.reshape(self.nx+1,self.ny,self.nz)
        
        hj = hi + (self.nx)*(self.ny+1)*self.nz
        ey = u[hi:hj]
        ey = ey.reshape(self.nx,self.ny+1,self.nz)
        
        ez = u[hj:]
        ez = ez.reshape(self.nx,self.ny,self.nz+1)
        
        return [ex,ey,ez]
    
    def pointSource(self, x,y,z):
        """ A routine to add a point source at the grid loc (x,y) """
        rhsz = np.zeros((self.nx,self.ny,self.nz+1),dtype='complex128') 
        rhsz[x,y,z] = 1.0
        rhsz = rhsz.flatten()
        rhsx = np.zeros((self.nx+1,self.ny,self.nz))
        rhsy = np.zeros((self.nx,self.ny+1,self.nz))

        self.rhs = np.concatenate((rhsx.flatten(), rhsy.flatten(), rhsz.flatten()))
    
    def planeWave(self):
        """ A routine to add a te planewave at angle as spec'd """
        thi = self.incAng 
        instep = 3+self.npml;
        # mdpt = nx/2; # should replace by div
        x = np.arange(1,1+self.nx,dtype='float64')*self.dx
        # The assumption is that the Ez and materials are co
        # located. Since epsilon(50) => in the half space, epsilon(51) =>
        # is not, the actual zero boundary must be between them, or on
        # that y boundary.
        # Mapping is y,x because of how matlab treats these. annoying.
        # (
        # print self.div+1
        Y,X = np.meshgrid(x-(self.div+1+0.5)*self.dx, x); # I do things backward
        
        Yyh, Xyh = np.meshgrid(np.append(0.0,x)+(self.dx/2)\
                               - (self.div+1+0.5)*self.dx,x);
        Yxh, Xxh = np.meshgrid(x-(self.div+1+0.5)*self.dx, \
                               np.append(0.0,x)+(self.dx/2));
                                 
        # matOut.savemat('grids', {'Y':Y, 'X':X, 'Yyh':Yyh, 'Xyh':Xyh, 'Yxh':Yxh, 'Xxh':Xxh})
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
        etaH = np.sqrt(self.muo/(self.eHS+(1j*self.sHS/self.w)));
        rTE = (ni*np.cos(thi) - nt*np.cos(tht))/(ni*np.cos(thi) + nt*np.cos(tht));
        tTE = (2*ni*np.cos(thi))/(ni*np.cos(thi) + nt*np.cos(tht));
          
        # print rTE
        # print tTE
    
        #   % [yy xx] = meshgrid(zeY, zeX);
        #   % Make a selector for the half space.
        ths = np.zeros([self.nx,self.ny],dtype='bool')
        ths[Y<0] = 1
        # ths = ths.astype(bool)
        #   size(ths);
    
        thsy = np.zeros([self.nx,self.nx+1],dtype='bool')
        thsy[Yyh<0] = 1
        # thsy = thsy.astype(bool)
    
        thsx = np.zeros([self.nx+1,self.nx],dtype='bool')
        thsx[Yxh<0] = 1
        #     thsx = thsx.astype(bool)
          
        # matOut.savemat('selctors', {'ths':ths, 'thsy':thsy, 'thsx':thsx})
    
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
          
        xl = instep-1; xr = self.nx-instep-1;
        yb = instep-1; yt = self.ny-instep-1;
        Jsrcz = np.zeros([self.nx,self.ny],dtype='complex128');
        Msrcx = np.zeros([self.nx,self.ny+1], dtype='complex128');
        Msrcy = np.zeros([self.nx+1,self.ny], dtype='complex128');
    
        Jsrcz[xl,yb:(yt+1)] =    Jsrcz[xl,yb:(yt+1)] + (1)*(Hyinc[xl,yb:(yt+1)]/self.dx);
        Jsrcz[xr,yb:(yt+1)] =    Jsrcz[xr,yb:(yt+1)] - (1)*(Hyinc[xr+1,yb:(yt+1)]/self.dx);
        Jsrcz[xl:(xr+1),yb] =    Jsrcz[xl:(xr+1),yb] - (1)*(Hxinc[xl:(xr+1),yb]/self.dy);
        Jsrcz[xl:(xr+1),yt] =    Jsrcz[xl:(xr+1),yt] + (1)*(Hxinc[xl:(xr+1),yt+1]/self.dy);
    
        Msrcx[xl:(xr+1),yb]   =  (1)*(Ezinc[xl:(xr+1),yb]/self.dy);
        Msrcx[xl:(xr+1),yt+1] = -(1)*(Ezinc[xl:(xr+1),yt]/self.dy);
      
        Msrcy[xl,   yb:(yt+1)] = -(1)*(Ezinc[xl,yb:(yt+1)]/self.dx);
        Msrcy[xr+1, yb:(yt+1)] =  (1)*(Ezinc[xr,yb:(yt+1)]/self.dx);
          

        pw = (1j*self.w*self.muo)*Jsrcz.flatten() - \
              sparse.kron(sparse.eye(self.nx,self.nx), self.d2/self.dx)*Msrcx.flatten() + \
              sparse.kron(self.d2/self.dx, sparse.eye(self.ny,self.ny))*Msrcy.flatten()
    
        self.rhs = pw.flatten();
        
        
    def getS(self):
        ''' return the coefficient necessary in the Md*P part to make things work '''
        return self.w*self.muo*1j
    
        