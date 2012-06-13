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
        
        self.N = 2*(self.nx+1)*self.nx + (self.nx+1)*(self.nx+1)
        self.sol = [np.zeros((self.nx,self.ny)), np.zeros((self.nx,self.ny))]

        for x in range(2):
            self.epsmap[x][:,:(div+1)] = self.eHS
            self.sigmap[x][:,:(div+1)] = self.sHS
    
    def parseFields(self, u):
        ex = u[xrange((self.nx+1)*(self.ny))]
        
        ex = ex.reshape(self.nx+1,self.ny)
        
            
            
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
        
    def getk(self):
        ''' method to return the diagonal for the self.nabla2 operator '''
        
        
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
        assert Mdx.shape == Mdy.shape
        
        
        
#        assert self.nRx*self.nRy == 
        self.Md = spt.vCat([Mdx, Mdy, sparse.coo_matrix((N,Mdx.shape[1]))]) #oh bugger, problems. 
        self.Md = self.Md.tocsc()
        self.Md = self.Md.T
#                              
                              
    def planeWave(self):
        ''' Make the plane wave illumination for the given frequency and inc angle '''
        instep = 3+self.npml;
        
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
  
        ni = 1
        nt = np.sqrt((self.eHS-self.sHS/(1j*self.w))*self.muo)/np.sqrt(self.epso*self.muo)
        #incident angle
        tht = np.arcsin(ni*np.sin(self.incAng)/nt);
#
#%   % create the coefficients to specify the space. 
#  kinc = -[sin(thi); cos(thi)]
#  ktx = -[sin(tht); cos(tht)];
#  kFS = 1i*sqrt(muo*epso*w^2);
#  etaF = sqrt(muo/epso);
#%   % [kFS kfs]
#  kHS = 1i*sqrt(muo*eHS*w^2 + 1i*w*muo*sHS);
#  etaH = sqrt(muo/(eHS+1i*sHS/w));
#%   rTE = (ni*cos(thi) - nt*cos(tht))/(ni*cos(thi) + nt*cos(tht));
#%   tTE = (2*ni*cos(thi))/(ni*cos(thi) + nt*cos(tht));
#
#rTM = (nt*cos(thi) - ni*cos(tht))/(ni*cos(tht) + nt*cos(thi));
#tTM = (2*ni*cos(thi))/(ni*cos(tht) + nt*cos(thi));
#
#%   % [yy xx] = meshgrid(zeY, zeX);
#%   % Make a selector for the half space.
#  ths = zeros(nx+1,ny+1);
#  res.sigHz = zeros(nx+1,nx+1);
#  ths(Yhh<=0) = 1;
#  ths = logical(ths);
#  size(ths);
#  res.sigHz(Yhh<=0) = sHS;
#  
#  thsy = zeros(nx,nx+1);
#  res.sigEy = zeros(nx,nx+1);
#  thsy(Yyh<=0) = 1;
#  thsy = logical(thsy);
#  res.sigEy(Yyh<=0) = sHS;
#  
#  thsx = zeros(nx+1,nx);
#  res.sigEx = zeros(nx+1,nx);
#  thsx(Yxh<=0) = 1;
#  thsx = logical(thsx);
#  res.sigEx(Yxh<=0) = sHS;
#  
#% % kHS
#% %  size(thsy)
#  Hzinc = zeros(nx+1,ny+1);
#  Hzinc(~ths) = (1/etaF)*exp(kFS*(Xhh(~ths)*kinc(1) + Yhh(~ths)*kinc(2))) + ...
#      rTM*(1/etaF)*exp(kFS*(Xhh(~ths)*kinc(1) - Yhh(~ths)*kinc(2)));
#  
#  Hzinc(ths) = tTM*(1/etaH)*exp(kHS*(Xhh(ths)*ktx(1) + Yhh(ths)*ktx(2)));
#  
#  Exinc = zeros(nx+1,nx);
#  Exinc(~thsx) = (kinc(2)*exp(kFS*(Xxh(~thsx)*kinc(1) + Yxh(~thsx)*kinc(2)))) + ...
#      rTM*(-kinc(2)*exp(kFS*(Xxh(~thsx)*kinc(1) - Yxh(~thsx)*kinc(2))));
#  
#  Exinc(thsx) = tTM*(ktx(2)*exp(kHS*(Xxh(thsx)*ktx(1) + Yxh(thsx)*ktx(2))));
#  Exinc = -Exinc;
#  
#  Eyinc = zeros(nx,nx+1);
#  Eyinc(~thsy) = (kinc(1)*exp(kFS*(Xyh(~thsy)*kinc(1) + ...
#                                   Yyh(~thsy)*kinc(2)))) + ...
#      rTM*(kinc(1)*exp(kFS*(Xyh(~thsy)*kinc(1) - ...
#                            Yyh(~thsy)*kinc(2))));
#  Eyinc(thsy) = tTM*(ktx(1)*exp(kHS*(Xyh(thsy)*ktx(1) + ...
#                                          Yyh(thsy)*ktx(2))));
#  
#  % Eyinc = -Eyinc;
#  
#  fldz.Hz = Hzinc;
#  fldz.Ex = Exinc;
#  fldz.Ey = Eyinc;
#  
#  fldz.sEx = reshape(hz2ex*Hzinc(:), nx+1,nx)./(-1i*w*epso+res.sigEx);
#  fldz.sEy = reshape(hz2ey*Hzinc(:), nx,nx+1)./(-1i*w*epso+res.sigEy);
#  fldz.sHz = reshape(ex2hz*Exinc(:) + ey2hz*Eyinc(:), nx+1,nx+1)/(1i*w*muo);
#  
#  % A = [spalloc(n,n,0) spalloc(n,n,0) hz2ex;
#  %    spalloc(n,n,0) spalloc(n,n,0) hz2ey;
#  %    ex2hz ey2hz spalloc(N,N,0)];
#  
#%   Hyinc = -Hyinc;
#%  %  HxincBB = (kron(d1/dx, speye(nx)) * Ezinc(:))/(1i*w*muo);
#%  %  HyincBB = (kron(speye(nx), -d1/dx) *Ezinc(:))/(1i*w*muo);
#%  %  Hxinc = reshape(Hxinc, nx,nx+1);
#  
#%  %  HxincBB = reshape(HxincBB,nx,nx+1);
#%  %  HyincBB = reshape(HyincBB, nx+1,nx);
#%  % % Debugging plots
#%  %  figure(18)
#%  %  subplot(211)
#%  %  plot(1:nx+1, real(Hxinc(50,:)), 1:nx+1, real(HxincBB(50,:)))
#%  %  subplot(212)
#%  %  plot(1:nx+1, imag(Hxinc(50,:)), 1:nx+1, imag(HxincBB(50,:)))
#  
#%  %  figure(19)
#%  %  subplot(211)
#%  %  plot(1:nx, real(Hyinc(50,:)), 1:nx, real(HyincBB(50,:)))
#%  %  subplot(212)
#%  %  plot(1:nx, imag(Hyinc(50,:)), 1:nx, imag(HyincBB(50,:)))
#  
#%  %    figure(2024)
#%  %  plot(real(Hxinc(50,:)))
#%  %  title('Hx at 50th row') - there is a discontinuity.
#  
#%  %  figure(2000);
#%  %  subplot(321)
#%  %  imagesc(real(Hxinc)')
#%  %  colorbar
#%  %  caxis([-4e-3 4e-3])
#%  %  cbgb=get(gca,'CLim')
#%  %  set(gca, 'YDir', 'normal')
#%  %  subplot(322)
#%  %  imagesc(real(HxincBB)')
#%  %  colorbar
#%  %  caxis(cbgb)
#  
#%  %  set(gca, 'YDir', 'normal')
#%  %    subplot(323)
#%  %  imagesc(real(Hyinc)')
#%  %  colorbar
#%  %  caxis([-3e-3 3e-3])
#%  %  cbgb=get(gca,'CLim');
#  
#%  %  set(gca, 'YDir', 'normal')
#%  %  subplot(324)
#%  %  imagesc(real(HyincBB)')
#%  %  colorbar
#%  %  caxis(cbgb)
#%  %  set(gca, 'YDir', 'normal')  
#%  %  subplot(325)
#%  %  imagesc(real(Ezinc)')
#%  %  colorbar
#%  %  set(gca, 'YDir', 'normal')
#%  %  subplot(326)
#%  %  imagesc(imag(Ezinc)')
#%  %  colorbar
#%  %  set(gca, 'YDir', 'normal')
#  
#%   xl = instep; xr = nx-instep;
#%   yb = instep; yt = ny-instep;
#
#%   Jsrcz = zeros(nx,ny);
#
#%   Msrcx = zeros(nx,ny+1);
#%   Msrcy = zeros(nx+1,ny);
#
#%   Jsrcz(xl,yb:yt) =    Jsrcz(xl,yb:yt) + (1)*(Hyinc(xl,yb:yt)/dx);
#%   Jsrcz(xr,yb:yt) =    Jsrcz(xr,yb:yt) - (1)*(Hyinc(xr+1,yb:yt)/dx);
#%   Jsrcz(xl:xr,yb) =    Jsrcz(xl:xr,yb) - (1)*(Hxinc(xl:xr,yb)/dy);
#%   Jsrcz(xl:xr,yt) =    Jsrcz(xl:xr,yt) + (1)*(Hxinc(xl:xr,yt+1)/dy);
#
#    
#%   Msrcx(xl:xr,yb)   =  (1)*(Ezinc(xl:xr,yb)/dy);
#%   Msrcx(xl:xr,yt+1) = -(1)*(Ezinc(xl:xr,yt)/dy);
#    
#%   Msrcy(xl,   yb:yt) = -(1)*(Ezinc(xl,yb:yt)/dx);
#%   Msrcy(xr+1, yb:yt) =  (1)*(Ezinc(xr,yb:yt)/dx);
#
#
#Msrcz = zeros(nx+1,nx+1);
#
#Jsrcx = zeros(nx+1,nx);
#Jsrcy = zeros(nx,nx+1);
#
#  instep = npml+3;
#  xl = instep; xr = nx-instep;
#  yb = instep; yt = ny-instep;
#  
#    % 5.48a
#  Jsrcx(xl+1:xr,yb) = Jsrcx(xl+1:xr,yb) + Hzinc(xl+1:xr,yb)/dx;
#  % 5.49a
#  Jsrcx(xl+1:xr,yt) = Jsrcx(xl+1:xr,yt) - Hzinc(xl+1:xr,yt+1)/dx;
#  
#  % 5.52a
#  Jsrcy(xl,yb+1:yt) = Jsrcy(xl,yb+1:yt) - Hzinc(xl,yb+1:yt)/dx;
#  % 5.53a
#  Jsrcy(xr,yb+1:yt) = Jsrcy(xr,yb+1:yt) + Hzinc(xr+1,yb+1:yt)/dx;
#  
#  % 5.54a
#  Msrcz(xl+1:xr,yb) = Msrcz(xl+1:xr,yb) - Exinc(xl+1:xr,yb)/dx;
#  
#  % 5.55a
#  Msrcz(xl+1:xr,yt+1) = Msrcz(xl+1:xr,yt+1) + Exinc(xl+1:xr,yt)/dx;
#    
#  % 5.58a
#  Msrcz(xl, yb+1:yt) = Msrcz(xl, yb+1:yt) + (1)*Eyinc(xl, yb+1:yt)/dx;
#  
#  % 5.59a
#  Msrcz(xr+1,yb+1:yt) = Msrcz(xr+1,yb+1:yt) - (1)*Eyinc(xr, yb+1:yt)/dx;
#
#  pw = [Jsrcx(:); Jsrcy(:); Msrcz(:)];

        
        