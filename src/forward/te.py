'''
Created on Jun 12, 2012

@author: dstrauss
'''

from model import fwd
import numpy as np
from scipy import sparse

class solver(fwd):
    ''' class to implement the transverse electric mode, rather - the case where we have Ez only ''' 
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