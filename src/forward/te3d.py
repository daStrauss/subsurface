'''
Created on Oct 10, 2012

@author: dstrauss
'''

from model import fwd
import numpy as np
from scipy import sparse
import sparseTools as spt
from scipy import io as spio

def speye(n):
    return sparse.eye(n,n)

class solver(fwd):
    ''' class to implement the transverse electric mode, rather - the case where we have Ez only ''' 
    def setspace(self, nx,ny,nz,dx, dy, dz):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.nz = nz # number z points
        self.N = (nx+1)*ny*nz + nx*(ny+1)*nz + nx*ny*(nz+1) # super size of space
        self.dx = dx # delta x
        self.dy = dy # delta y
        self.dz = dz # delta z

        self.npml = min(10,round(nx/10))
    
    def makeGradOper(self):
        ''' routine to make a big matrix for TE problems ex,ey,ez all incorporated,
        based on the petsc_cpx routine.''' 
        # a quick hint: pd2 == pdo
        pd1 = self.ph*self.d1 # maps full to half grid
        pd2 = self.po*self.d2 # maps half to full
        
        AA = sparse.kron(speye(self.nx+1), sparse.kron(pd2,speye(self.nz)))*\
             sparse.kron(speye(self.nx+1), sparse.kron(pd1,speye(self.nz))) + \
             sparse.kron(speye(self.nx+1), sparse.kron(speye(self.ny),pd2)) * \
             sparse.kron(speye(self.nx+1), sparse.kron(speye(self.ny),pd1)) 
# chngd
        AB = sparse.kron(speye(self.nx+1), sparse.kron(pd2,speye(self.nz)))*\
             sparse.kron(pd1, sparse.kron(speye(self.ny+1),speye(self.nz)))
            
        AC = sparse.kron(speye(self.nx+1),sparse.kron(speye(self.ny),pd2))*\
            sparse.kron(pd1,sparse.kron(speye(self.ny),speye(self.nz+1)))
# chngd
        BA = sparse.kron(pd2,sparse.kron(speye(self.ny+1),speye(self.nz)))*\
             sparse.kron(speye(self.nx+1),sparse.kron(pd1,speye(self.nz)))
# chngd
        BB = sparse.kron(pd2,sparse.kron(speye(self.ny+1),speye(self.nz)))*\
             sparse.kron(pd1,sparse.kron(speye(self.ny+1),speye(self.nz))) + \
             sparse.kron(speye(self.nx),sparse.kron(speye(self.ny+1),pd2))*\
             sparse.kron(speye(self.nx),sparse.kron(speye(self.ny+1),pd1))
# chngd
        BC = sparse.kron(speye(self.nx),sparse.kron(speye(self.ny+1),pd2))*\
             sparse.kron(speye(self.nx),sparse.kron(pd1,speye(self.nz+1)))
# chngd
        CA = sparse.kron(pd2,sparse.kron(speye(self.ny),speye(self.nz+1)))*\
             sparse.kron(speye(self.nx+1),sparse.kron(speye(self.ny),pd1))
# chngd
        CB = sparse.kron(speye(self.nx),sparse.kron(pd2,speye(self.nz+1)))*\
            sparse.kron(speye(self.nx),sparse.kron(speye(self.ny+1),pd1))
# chngd
        CC = sparse.kron(speye(self.nx),sparse.kron(pd2, speye(self.nz+1)))*\
             sparse.kron(speye(self.nx),sparse.kron(pd1, speye(self.nz+1))) + \
             sparse.kron(pd2,sparse.kron(speye(self.ny),speye(self.nz+1)))*\
             sparse.kron(pd1,sparse.kron(speye(self.ny),speye(self.nz+1)))
# chngd
            
            # legacy - matlab ordering 
#        AA = sparse.kron(speye(self.nz),sparse.kron(pd2,speye(self.nx+1)))*\
#            sparse.kron(speye(self.nz),sparse.kron(pd1,speye(self.nx+1))) + \
#            sparse.kron(pd2,speye((self.nx+1)*self.ny))*sparse.kron(pd1,speye((self.nx+1)*self.ny))
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
        sigX[:,:(div),:] = self.sHS
        
        sigY = np.zeros((self.nx,self.ny+1,self.nz))
        sigY[:,:(div+1),:] = self.sHS
        
        sigZ = np.zeros((self.nx,self.ny,self.nz+1))
        sigZ[:,:(div),:] = self.sHS
        
        self.sigX = sigX
        self.sigY = sigY
        self.sigZ = sigZ
        
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
        '''Creates an n-grid mesh across the surface for the 3D case '''
        
        self.nSen = nSensors*nSensors
        '''First find the appropriate 10 indexes within the PML & illumination region '''
        indx = np.round(np.linspace(self.npml+5,self.nx-self.npml-5, nSensors)).astype(int)-1;
        indx = np.unique(indx)
        # print (indx + 1)
        ''' make the exact X operator using strides '''
        xl,zl = np.meshgrid(indx+1,indx)
        Mx = sparse.dok_matrix((self.nSen,(self.nx+1)*self.ny*self.nz))
        
        for ix,loc in enumerate(zip(xl.flatten(),zl.flatten())):
            pts = loc[0]*self.ny*self.nz + self.div*self.nz + loc[1]
            Mx[ix,pts] = 1.0
        
        xl,zl = np.meshgrid(indx,indx)
        My = sparse.dok_matrix((self.nSen,self.nx*(self.ny+1)*self.nz))
        
        for ix,loc in enumerate(zip(xl.flatten(),zl.flatten())):
            pts = loc[0]*(self.ny+1)*self.nz + (self.div+1)*self.nz + loc[1]
            My[ix,pts] = 1.0
            
            
        '''make the exact Z operator using strides '''
        xl,zl = np.meshgrid(indx,indx+1)
        Mz = sparse.dok_matrix((self.nSen,self.nx*self.ny*(self.nz+1)))
        
        for ix,loc in enumerate(zip(xl.flatten(),zl.flatten())): 
            pts = loc[0]*self.ny*(self.nz+1) + self.div*(self.nz+1) + loc[1]
            Mz[ix,pts] = 1.0        
        
        ''' smush together in block diagonal format '''
        self.Ms = sparse.block_diag((Mx,My,Mz),'csr')
        
    def setCTRX(self):
        ''' create some operators to map back and forth between the x space and the u '''
        self.p2x = sparse.eye(self.nRx*self.nRy*self.nRz,self.nRx*self.nRy*self.nRz)
        self.p2x = sparse.vstack((self.p2x,self.p2x,self.p2x))
        
        print 'p2x size  ' + repr(self.p2x.shape)
        if not(hasattr(self, 'x2u')):
            print 'Oh youre hosed now!'
        
#        # print self.p2x.shape
#        self.x2u = self.Md.T
#        # print self.x2u.shape
    
    def getXSize(self):
        ''' return the proper size of X so that the optimization routine can work its magic '''
        return 3*self.nRx*self.nRy*self.nRz
    def getPSize(self):
        ''' remember, it changes from solver to solver '''
        return self.nRx*self.nRy*self.nRz
    
    def parseP(self,P):
        ''' the other two solvers are simple, this one is different '''
        return P.reshape(self.nRx,self.nRy,self.nRz)
        
        
    def setMd(self, xrng, yrng, zrng):
        '''Tell me the xrange,yrange, and zrange and Ill
        1) specify nRx,nRy, and nRz
        2) produce a matrix that achieves a 1:1 sampling, self.Md '''
        
        '''set the right dimensions'''
        self.nRx = xrng[1]-xrng[0]
        self.nRy = yrng[1]-yrng[0]
        self.nRz = zrng[1]-zrng[0]
        
        nR = self.nRx*self.nRy*self.nRz
        ''' ok have to use spans:
        loc = i*J*K + j*K + k for row-major ordering '''
        ''' populate the locations in the X grid'''
        #sX = sparse.dok_matrix((self.nx+1,self.ny,self.nz),dtype='bool')
        #sX[xrng[0]+1:xrng[1]+1,yrng[0]:yrng[1],zrng[0]:zrng[1]] = True
        ''' make it an operator '''
        ''' nested for should give reshape-able vectors '''
        cnt = 0
        Mx = sparse.dok_matrix((nR,(self.nx+1)*self.ny*self.nz))
        for x in xrange(xrng[0]+1,xrng[1]+1):
            for y in xrange(yrng[0],yrng[1]):
                for z in xrange(zrng[0],zrng[1]):
                    pts = x*self.ny*self.nz + y*self.nz + z
                    Mx[cnt,pts] = 1.0
                    cnt += 1
        
        '''populate the locations in the Y grid'''
        My = sparse.dok_matrix((nR,self.nx*(self.ny+1)*self.nz))
        cnt = 0
        for x in xrange(xrng[0],xrng[1]):
            for y in xrange(yrng[0]+1,yrng[1]+1):
                for z in xrange(zrng[0],zrng[1]):
                    pts = x*(self.ny+1)*self.nz + y*self.nz + z
                    My[cnt,pts] = 1.0
                    cnt += 1
        
        
        '''populate the locations in the Z grid'''
        Mz = sparse.dok_matrix((nR,self.nx*self.ny*(self.nz+1)))
        cnt = 0
        for x in xrange(xrng[0],xrng[1]):
            for y in xrange(yrng[0],yrng[1]):
                for z in xrange(zrng[0]+1,zrng[1]+1):
                    pts = x*(self.ny)*(self.nz+1) + y*(self.nz+1) + z
                    Mz[cnt,pts] = 1.0
                    cnt += 1
        
        ''' put them all together in a block matrix '''    
        self.Md = spt.vCat([Mx.T,My.T,Mz.T]).T
        print 'Md shape ' + repr(self.Md.shape)
        
        self.x2u = sparse.block_diag((Mx,My,Mz), 'csc').T
        print 'x2u shape ' + repr(self.x2u.shape)
        

        
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
        # rhsz = rhsz.flatten()
        rhsx = np.zeros((self.nx+1,self.ny,self.nz))
        rhsy = np.zeros((self.nx,self.ny+1,self.nz))
        
        #self.rhs = np.zeros(self.N,dtype='complex128')
        #self.rhs[23614-1] = 1.0
        self.rhs = np.concatenate((rhsx.flatten(), rhsy.flatten(), rhsz.flatten()))
    
    def planeWave(self):
        ''' populates self.rhs with a planewave with incident conditions:
        self.incAng, self.azAng '''
        
        ''' incident angles, local copy '''
        thi = self.incAng
        phi = self.azAng 
        
        ''' how far in from the PML should we go? -- 2 grid points should be enough '''
        instep = 2+self.npml;
        # mdpt = nx/2; # should replace by div
        x = np.arange(1,1+self.nx,dtype='float64')*self.dx
        y = np.arange(1,1+self.ny,dtype='float64')*self.dy
        z = np.arange(1,1+self.nz,dtype='float64')*self.dz
        # The assumption is that the Ez and materials are co
        # located. Since epsilon(50) => in the half space, epsilon(51) =>
        # is not, the actual zero boundary must be between them, or on
        # that y boundary.
        # Mapping is y,x because of how matlab treats these. annoying.
        # (
        # print self.div+1
        Yz,Xz,Zz = np.meshgrid(y-(self.div+0.5)*self.dy,\
                               x,\
                               (np.append(0.0,z) + self.dz/2.0))
        Yx,Xx,Zx = np.meshgrid(y-(self.div+0.5)*self.dy,\
                               (np.append(0.0,x) + self.dx/2.0),\
                               z)
        Yy,Xy,Zy = np.meshgrid((np.append(0.0,y) + self.dy/2.0)-(self.div+0.5)*self.dy,\
                               x,\
                               z)
        
        Yxh,Xxh,Zxh = np.meshgrid(np.append(0.0,y)+self.dy/2.0 - (self.div+0.5)*self.dy, \
                                  x,\
                                  np.append(0.0,z)+self.dz/2.0)
        
        Yyh,Xyh,Zyh = np.meshgrid(y-(self.div+0.5)*self.dy,\
                                  np.append(0.0,x)+self.dx/2.0,\
                                  np.append(0.0,z)+self.dz/2.0)
        
        Yzh,Xzh,Zzh = np.meshgrid(np.append(0.0,y)+self.dy/2.0-(self.div+0.5)*self.dy,\
                                  np.append(0.0,x)+self.dx/2.0,\
                                  z)
        
                                 
        # matOut.savemat('grids', {'Y':Y, 'X':X, 'Yyh':Yyh, 'Xyh':Xyh, 'Yxh':Yxh, 'Xxh':Xxh})
        ni = 1;
        nt = np.sqrt((self.eHS-self.sHS/(1j*self.w))*self.muo)/np.sqrt(self.epso*self.muo);
    
        # transmitted angle angle
        # thi = 45*np.pi/180 # taken as input argument.
        tht = np.arcsin(ni*np.sin(thi)/nt);
    
        # create the coefficients to specify the space. 
        kinc = -1*np.array([np.sin(thi)*np.cos(phi), np.cos(thi), -np.sin(phi)*np.sin(thi)])
        ktx =  -1*np.array([np.sin(tht)*np.cos(phi), np.cos(tht), -np.sin(phi)*np.sin(thi)])
        kFS = 1j*self.kfree;
        kHS = 1j*self.kHS;
        
        etaF = np.sqrt(self.muo/self.epso);
        etaH = np.sqrt(self.muo/(self.eHS+(1j*self.sHS/self.w)));
        
        rTE = (ni*np.cos(thi) - nt*np.cos(tht))/(ni*np.cos(thi) + nt*np.cos(tht));
        tTE = (2*ni*np.cos(thi))/(ni*np.cos(thi) + nt*np.cos(tht));
          
        ''' make selectors for the subspace and so on'''
        thsz = np.zeros((self.nx,self.ny,self.nz+1),dtype='bool')
        thsz[Yz<0] = 1

        thsx = np.zeros((self.nx+1,self.ny,self.nz),dtype='bool')
        thsx[Yx<0] = 1
        # thsy = thsy.astype(bool)
    
        thsy = np.zeros((self.nx,self.ny+1,self.nz),dtype='bool')
        thsy[Yy<0] = 1
        
        ''' make selectors for the halfgrid subspace '''
        
        thsxh = np.zeros((self.nx,self.ny+1,self.nz+1),dtype='bool')
        thsxh[Yxh<0] = 1
        
        thsyh = np.zeros((self.nx+1,self.ny,self.nz+1),dtype='bool')
        thsyh[Yyh<0] = 1
        
        thszh = np.zeros((self.nx+1,self.ny+1,self.nz),dtype='bool')
        thszh[Yzh<0] = 1
        
        
        ''' in the matlab version, I stick some eHS in here at this point - not sure if I need that now '''
        
        Ezinc = np.cos(phi)*te_ezf(Xz,Yz,Zz,thsz,kinc,ktx,rTE,tTE,kFS,kHS)
        Exinc = np.sin(phi)*te_ezf(Xx,Yx,Zx,thsx,kinc,ktx,rTE,tTE,kFS,kHS)
        Eyinc = np.zeros((self.nx,self.ny+1,self.nz),dtype='complex128')
        
        
        
        Hxinc = te_hf(Xxh,Yxh,Zxh,thsxh,kinc,ktx,rTE,tTE,kFS,kHS)
        Hxinc[~thsxh] = Hxinc[~thsxh]*(1.0/etaF)*np.cos(phi)*kinc[1]
        Hxinc[thsxh] = Hxinc[thsxh]*(1.0/etaH)*np.cos(phi)*ktx[1]
        
        Hyinc = te_ezf(Xyh,Yyh,Zyh,thsyh,kinc,ktx,rTE,tTE,kFS,kHS)
        Hyinc[~thsyh] = Hyinc[~thsyh]*(1.0/etaF)*(kinc[2]*np.sin(phi) - kinc[0]*np.cos(phi))
        Hyinc[thsyh] = Hyinc[thsyh]*(1.0/etaH)*(ktx[2]*np.sin(phi) - ktx[0]*np.cos(phi))
        
        Hzinc = te_ezf(Xzh,Yzh,Zzh,thszh,kinc,ktx,rTE,tTE,kFS,kHS)
        Hzinc[~thszh] = Hzinc[~thszh]*(1.0/etaF)*(-np.sin(phi)*kinc[1])
        Hzinc[thszh] = Hzinc[thszh]*(1.0/etaH)*(-np.sin(phi)*kinc[1])
        
        # spio.savemat('intrn', {'hxinc':Hxinc,'hyinc':Hyinc,'hzinc':Hzinc})
          
        xl = instep-1; xr = self.nx-instep-1;
        yb = instep-1; yt = self.ny-instep-1;
        zb = instep-1; zt = self.nz-instep-1;
        
        Jsrcx = np.zeros([self.nx+1,self.ny,self.nz],dtype='complex128')
        Jsrcy = np.zeros([self.nx,self.ny+1,self.nz],dtype='complex128')
        Jsrcz = np.zeros([self.nx,self.ny,self.nz+1],dtype='complex128')
        
        Msrcx = np.zeros([self.nx,self.ny+1,self.nz+1], dtype='complex128')
        Msrcy = np.zeros([self.nx+1,self.ny,self.nz+1], dtype='complex128')
        Msrcz = np.zeros([self.nx+1,self.ny+1,self.nz], dtype='complex128')
        
         # 5.48a
        Jsrcx[xl+1:(1+xr),yb,zb:(1+zt)] = Jsrcx[xl+1:(1+xr),yb,zb:(1+zt)] + Hzinc[xl+1:(1+xr),yb,zb:(1+zt)]/self.dx
        
        # 5.49a
        Jsrcx[xl+1:(1+xr),yt,zb:(1+zt)] = Jsrcx[xl+1:(1+xr),yt,zb:(1+zt)] - Hzinc[xl+1:(1+xr),yt+1,zb:(1+zt)]/self.dx
  
        # 5.48b
        Jsrcz[xl:(1+xr),yb,zb+1:(1+zt)] = Jsrcz[xl:(1+xr),yb,zb+1:(1+zt)] - Hxinc[xl:(1+xr),yb,zb+1:(1+zt)]/self.dx
        
        # 5.49b
        Jsrcz[xl:(1+xr),yt,zb+1:(1+zt)] = Jsrcz[xl:(1+xr),yt,zb+1:(1+zt)] + Hxinc[xl:(1+xr),yt+1,zb+1:(1+zt)]/self.dx
        
        # 5.50a
        Jsrcx[xl+1:(1+xr),yb:(1+yt),zb] = Jsrcx[xl+1:(1+xr),yb:(1+yt),zb] - Hyinc[xl+1:(1+xr),yb:(1+yt),zb]/self.dx
        
        # 5.50b
        Jsrcy[xl:(1+xr),yb+1:(1+yt),zb] = Jsrcy[xl:(1+xr),yb+1:(1+yt),zb] + Hxinc[xl:(1+xr),yb+1:(1+yt),zb]/self.dx
  
        # 5.51a
        Jsrcx[xl+1:(1+xr),yb:(1+yt),zt] = Jsrcx[xl+1:(1+xr),yb:(1+yt),zt] + Hyinc[xl+1:(1+xr),yb:(1+yt),zt+1]/self.dx
        
        # 5.51b
        Jsrcy[xl:(1+xr),yb+1:(1+yt),zt] = Jsrcy[xl:(1+xr),yb+1:(1+yt),zt] - Hxinc[xl:(1+xr),yb+1:(1+yt),zt+1]/self.dx
  
        # 5.52a
        Jsrcy[xl,yb+1:(1+yt),zb:(1+zt)] = Jsrcy[xl,yb+1:(1+yt),zb:(1+zt)] - Hzinc[xl,yb+1:(1+yt),zb:(1+zt)]/self.dx
        
        # 5.52b
        Jsrcz[xl,yb:(1+yt),zb+1:(1+zt)] = Jsrcz[xl,yb:(1+yt),zb+1:(1+zt)] + Hyinc[xl,yb:(1+yt),zb+1:(1+zt)]/self.dx
        
        # 5.53a
        Jsrcy[xr,yb+1:(1+yt),zb:(1+zt)] = Jsrcy[xr,yb+1:(1+yt),zb:(1+zt)] + Hzinc[xr+1,yb+1:(1+yt),zb:(1+zt)]/self.dx
  
        # 5.53b 
        Jsrcz[xr,yb:(1+yt),zb+1:(1+zt)] = Jsrcz[xr,yb:(1+yt),zb+1:(1+zt)] - Hyinc[xr+1,yb:(1+yt),zb+1:(1+zt)]/self.dx
        
        # 5.54a
        Msrcz[xl+1:(1+xr),yb,zb:(1+zt)] = Msrcz[xl+1:(1+xr),yb,zb:(1+zt)] - Exinc[xl+1:(1+xr),yb,zb:(1+zt)]/self.dx
  
        # 5.54b
        Msrcx[xl:(1+xr),yb,zb+1:(1+zt)] = Msrcx[xl:(1+xr),yb,zb+1:(1+zt)] + Ezinc[xl:(1+xr),yb,zb+1:(1+zt)]/self.dx
  
        # 5.55a
        Msrcz[xl+1:(1+xr),yt+1,zb+1:(1+zt)] = Msrcz[xl+1:(1+xr),yt+1,zb+1:(1+zt)] + Exinc[xl+1:(1+xr),yt, zb+1:(1+zt)]/self.dx
  
        # 5.55b
        Msrcx[xl:(1+xr),yt+1,zb+1:(1+zt)] = Msrcx[xl:(1+xr),yt+1,zb+1:(1+zt)] - Ezinc[xl:(1+xr),yt,zb+1:(1+zt)]/self.dx
  
        # 5.56a
        Msrcy[xl+1:(1+xr),yb:(1+yt),zb] = Msrcy[xl+1:(1+xr),yb:(1+yt),zb] + Exinc[xl+1:(1+xr),yb:(1+yt),zb]/self.dx
  
        # 5.56b
        Msrcx[xl:(1+xr),yb+1:(1+yt),zb] = Msrcx[xl:(1+xr),yb+1:(1+yt),zb] - Eyinc[xl:(1+xr),yb+1:(1+yt),zb]/self.dx
  
        # 5.57a
        Msrcy[xl+1:(1+xr),yb:(1+yt),zt+1] = Msrcy[xl+1:(1+xr),yb:(1+yt),zt+1] - Exinc[xl+1:(1+xr),yb:(1+yt),zt]/self.dx
  
        # 5.57b
        Msrcx[xl:(1+xr),yb+1:(1+yt),zt+1] = Msrcx[xl:(1+xr),yb+1:(1+yt),zt+1] + Eyinc[xl:(1+xr),yb+1:(1+yt),zt]/self.dx
  
        # 5.58a
        Msrcz[xl, yb+1:(1+yt),zb:(1+zt)] = Msrcz[xl, yb+1:(1+yt), zb:(1+zt)] + Eyinc[xl, yb+1:(1+yt), zb:(1+zt)]/self.dx
  
        # 5.58b
        Msrcy[xl, yb:(1+yt), zb+1:(1+zt)] = Msrcy[xl,yb:(1+yt),zb+1:(1+zt)] - Ezinc[xl,yb:(1+yt),zb+1:(1+zt)]/self.dx
  
        # 5.59a
        Msrcz[xr+1,yb+1:(1+yt),zb:(1+zt)] = Msrcz[xr+1,yb+1:(1+yt),zb:(1+zt)] - Eyinc[xr, yb+1:(1+yt),zb:(1+zt)]/self.dx
        
        # 5.59b
        Msrcy[xr+1,yb:(1+yt),zb+1:(1+zt)] =  Msrcy[xr+1,yb:(1+yt),zb+1:(1+zt)] + Ezinc[xr,yb:(1+yt),zb+1:(1+zt)]/self.dx
        # spio.savemat('intrn', {'hxinc':Hxinc,'hyinc':Hyinc,'hzinc':Hzinc})
          
        
        spio.savemat('intrn', {'sigx':self.sigX, 'sigy':self.sigY, 'sigz':self.sigZ,\
                               'exinc':Exinc,'eyinc':Eyinc,'ezinc':Ezinc,\
                               'hxinc':Hxinc,'hyinc':Hyinc,'hzinc':Hzinc,\
                               'jsrcx':Jsrcx,'jsrcy':Jsrcy,'jsrcz':Jsrcz, \
                               'msrcx':Msrcx,'msrcy':Msrcy,'msrcz':Msrcz})
          
        
        J = np.concatenate((Jsrcx.flatten(),Jsrcy.flatten(),Jsrcz.flatten()))
        M = np.concatenate((Msrcx.flatten(),Msrcy.flatten(),Msrcz.flatten()))
        
        # ok. have to make some transformations in order to make this happen.
        pd1 = self.ph*self.d1 # maps full to half grid
        pd2 = self.po*self.d2
        
        '''make a curl operator from M -> J '''
        nex = (self.nx+1)*self.ny*self.nz
        ney = self.nx*(self.ny+1)*self.nz
        nez = self.nx*self.ny*(self.nz+1)
        
        nhx = self.nx*(self.ny+1)*(self.nz+1)
        nhy = (self.nx+1)*self.ny*(self.nz+1)
        nhz = (self.nx+1)*(self.ny+1)*self.nz
                
        AA = sparse.coo_matrix((nex,nhx))
        
        AB = -sparse.kron(speye(self.nx+1),sparse.kron(speye(self.ny),pd2))
        
        AC = sparse.kron(speye(self.nx+1),sparse.kron(pd2,speye(self.nz)))
        
        BA = sparse.kron(speye(self.nx+1),sparse.kron(speye(self.ny),pd2))
        
        BB = sparse.coo_matrix((ney,nhy))
        
        BC = -sparse.kron(pd2,sparse.kron(speye(self.ny+1),speye(self.nz)))
        
        CA = -sparse.kron(speye(self.nx),sparse.kron(pd2,speye(self.nz+1)))
        
        CB = sparse.kron(pd2,sparse.kron(speye(self.ny),speye(self.nz+1)))
        
        CC = sparse.coo_matrix((nez,nhz))
        
        srB = spt.vCat([spt.hCat([AA,AB,AC]), spt.hCat([BA,BB,BC]), spt.hCat([CA,CB,CC])])
        
        self.rhs = -(1j*self.muo*self.w*J + srB*M)
        
    
    def getS(self):
        ''' return the coefficient necessary in the Md*P part to make things work '''
        return self.w*self.muo*1j
    
def te_ezf(X,Y,Z,xi,kinc,ktx,rTE,tTE,kFS,kHS):
    m,n,p = X.shape
    
    ez = np.zeros((m,n,p),dtype='complex128')
    ez[~xi] = np.exp(kFS*(X[~xi]*kinc[0]  + \
                          Y[~xi]*kinc[1]  + \
                          Z[~xi]*kinc[2])) + \
              rTE*np.exp(kFS*(X[~xi]*kinc[0] - \
                              Y[~xi]*kinc[1] + \
                              Z[~xi]*kinc[2]))
    ez[xi] = tTE*np.exp(kHS*(X[xi]*ktx[0]  + \
                          Y[xi]*ktx[1]  + \
                          Z[xi]*ktx[2]))
    return ez

    
def te_hf(X,Y,Z,xi,kinc,ktx,rTE,tTE,kFS,kHS):
    ''' right the only difference is in the sign of the reflected wave '''
    m,n,p = X.shape
    
    ez = np.zeros((m,n,p),dtype='complex128')
    ez[~xi] = np.exp(kFS*(X[~xi]*kinc[0]  + \
                          Y[~xi]*kinc[1]  + \
                          Z[~xi]*kinc[2])) - \
              rTE*np.exp(kFS*(X[~xi]*kinc[0] - \
                              Y[~xi]*kinc[1] + \
                              Z[~xi]*kinc[2]))
    ez[xi] = tTE*np.exp(kHS*(X[xi]*ktx[0]  + \
                          Y[xi]*ktx[1]  + \
                          Z[xi]*ktx[2]))
    return ez

    