'''
Created on Jun 12, 2012

@author: dstrauss
'''

from scipy import sparse
from scipy.sparse import linalg as lin
import scipy.special as spec
import pickle
import scipy.io as spio
import numpy as np

class pmlList(object):
    freq = np.ndarray(0)
    dx = np.ndarray(0)
    pmc = np.ndarray(0)

class fwd(object):
    epso = 8.854e-12
    muo = 4*np.pi*1e-7
    c = np.sqrt(1/(muo*epso))
    
    def __init__(self, freq, incAng):   
        '''initialize object ''' 
        self.f = freq
        self.w = 2*np.pi*freq
        self.l = self.c/self.f
        self.incAng = incAng
        
    def setspace(self, nx,ny,dx,dy):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.N = nx*ny # super size of space
        self.dx = dx # delta x
        self.dy = dy # delta y

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
        self.sol = [np.zeros((self.nx,self.ny)), np.zeros((self.nx,self.ny))]

        for x in range(2):
            self.epsmap[x][:,:(div+1)] = self.eHS
            self.sigmap[x][:,:(div+1)] = self.sHS
        
    def getk(self, ind):
        """ This routine assembles a diagonal matrix with the materials indexed by ind
        """
        kl = (self.muo*self.epsmap[ind]*(self.w**2) + 1j*self.w*self.muo*self.sigmap[ind]   )
        
        return sparse.spdiags(kl.flatten(), 0, self.N, self.N)

    
    def getPMLparm(self):
        ''' Returns the proper PML parameter either from a library if it has already been calculated,
        or it is able to explicitly calculate it otherwise
        '''
        try:
            lkt = pickle.load(open('pmlLib.p', 'rb'))
        except:
            lkt = pmlList()
            
        if np.any(lkt.freq==self.f):
            pmc = lkt.pmc[(lkt.freq==self.f) & (lkt.dx == self.dx)]
            # print pmc
            return pmc
        else:
            pmc = findBestAng(self.f, self.dx)
            lkt.freq = np.append(lkt.freq, self.f)
            lkt.pmc = np.append(lkt.pmc, pmc)
            lkt.dx = np.append(lkt.pmc, self.dx)
            pickle.dump(lkt, open('pmlLib.p', 'wb'))
            # print pmc
            return pmc
            
    def setOperators(self):
        ''' An wraper for makeFinteDiferences that will select the right PML parameters
        '''
        # self.makeFD(1886.7678)
        self.makeFD(self.getPMLparm())    
            
        
    def makeFD(self, pmc):
        """A routine that will define the following operators:
        A, ox, oy"""
        # h = (self.dx)*np.mat(np.ones((self.nx+1,1)));
        hi = np.zeros(self.nx+1);
    
        hi[:self.npml] = (pmc)*np.linspace(1,0,self.npml);
        hi[(self.nx-self.npml+1):(self.nx+1)] = \
                  (pmc)*np.linspace(0,1,self.npml);

        h = self.dx*np.ones(self.nx+1,dtype='complex128') + 1j*hi
        # h = self.dx*np.ones(self.nx+1,dtype='complex128')
        # matOut.savemat('Hout', {'h':h})
        
        opr = np.ones((2,self.nx+1))
        opr[0,:] = -1
        
        self.d1 = sparse.spdiags(opr, [-1, 0], self.nx+1, self.nx);
        self.d2 = sparse.spdiags(opr, [0,  1], self.nx, self.nx+1);

        # averaging operator for getting to the in-between grid points
        avg_n_c = sparse.spdiags(0.5*np.ones((2,self.nx+1)), [0, 1], self.nx,self.nx+1);

        #creating the half and non-half space stretchers.
        self.ph = sparse.spdiags(1/h, 0, self.nx+1,self.nx+1);
        self.po = sparse.spdiags(1/(avg_n_c*h), 0, self.nx,self.nx);
        
    def makeGradOper(self):
        ''' going to separate out the grad operator making so that I can override it in the TM case '''
        oper1d = (self.po*self.d2)*(self.ph*self.d1);
        # ok. building nabla2 for the TE only
        self.nabla2 = sparse.kron(oper1d,sparse.eye(self.nx,self.nx)) \
                 + sparse.kron(sparse.eye(self.nx,self.nx), oper1d);
                 
        
    def point_source(self, x,y):
        """ A routine to add a point source at the grid loc (x,y) """
        self.rhs.reshape(self.nx,self.ny)[x,y] = -1.0



    def fwd_solve(self, ind):
        """ Does the clean solve, prints a figure """
        # self.sol[ind] = self.rhs.copy();
        # b = self.rhs.copy().flatten();
        self.sol[ind] = lin.spsolve(sparse.csc_matrix(self.nabla2+self.getk(ind)),\
                                    self.rhs.flatten())
        # Q = lin.factorized(sparse.csc_matrix(self.nabla2 + self.getk(ind)))
        
        # self.sol[ind] = Q(self.rhs.flatten())
        # umfpack.linsolv((self.nabla2 + self.getk(ind)), self.sol[ind])
        # self.sol[ind] = np.array(self.sol[ind])
        # self.sol[ind] = self.sol[ind].reshape(self.nx,self.ny)
        
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
        
        
    def setMd(self, xrng, yrng):
        '''Tell me the xrange and the yrange and I'll make selector
        '''
        oprx = np.zeros((self.nx,self.ny),dtype='bool')
        oprx[xrng[0]:xrng[1],yrng[0]:yrng[1]] = 1
        self.nRx = xrng[1]-xrng[0]
        self.nRy = yrng[1]-yrng[0]
        
        idx = np.arange(self.N)
        oprx = oprx.flatten()
        idx = idx[oprx]
        self.Md = sparse.dok_matrix((self.N,idx.size))
        
        for i in range(idx.size):
            self.Md[idx[i],i]=1.0
        self.Md = self.Md.tocsc()
        self.Md = self.Md.T
        
    def plotSol(self,ind):
        pass

    def writeOut(self,ind):
        D = {'f':self.f, 'angle':self.incAng, 'sigMat':self.simap[ind], 'fld':self.sol[ind]}
        spio.savemat('maxPrint' + repr(self.rank), D)
                     
    
def findBestAng(freq,dx):
    '''for a particular frequency, find the best PML complex angle 
    '''
    # set some internal parameters
    nx = 99
    ny = nx
    # dx = 5.0
    dy = dx
    
    # this means free space
    eHS = 1.0
    sHS = 0
    
    epso = 8.854e-12
    muo = 4.0*np.pi*1e-7
    
    x = np.array(range(nx))*dx - dx*(nx/2)
    Y,X = np.meshgrid(x, x);
    dist = np.sqrt(Y**2+X**2) 
    
    c = 1/np.sqrt(muo*epso)
    
    k = (2*np.pi)/(c/freq)
    bbx = dx*dx*(1j/4)*spec.hankel1(0,k*dist)
    mdpt = nx/2
    # the center contains a singularity, but we'll  be forgiving and finite
    bbx[mdpt,mdpt] = 0.25*bbx[mdpt+1,mdpt] + 0.25*bbx[mdpt-1,mdpt] + 0.25*bbx[mdpt,mdpt+1] + 0.25*bbx[mdpt,mdpt-1]
    
    # sze = bbx.shape
    # print sze
    mask = np.ones(bbx.shape)
    mask[:15,:] =0
    mask[(nx-15):(nx),:] = 0
    mask[:,:15]=0
    mask[:,(nx-15):(nx)] = 0
    
#    print mask[:,98]
#    
#    plt.figure(388)
#    plt.imshow(mask)
#    plt.show()
#       
    lo = 0
    hi = 5
    
    # do two loops over different ranges to get a more precise estimate
    for tune in range(2):
        localError = np.zeros(50)
        angChoice = np.logspace(lo,hi,50)
    
        for i in range(50):    
            #     create a single object
            bce = fwd(freq)
            bce.setspace(nx,ny,dx,dy)
            bce.setmats(eHS, sHS, nx/2)
    
            # pmc = 22.229964825261945
            bce.makeFD(angChoice[i])
            bce.makeGradOper()
            bce.point_source(49, 49)
            bce.fwd_solve(0)
    
            localError[i] = np.linalg.norm((bce.sol[0] - bbx)*mask,'fro')
            # print 'ang = ' + repr(angChoice[i]) + ' local error ' + repr(localError[i])
    
        minIdx = np.argmin(localError)
        
        lo = np.log10(angChoice[max(0,minIdx-1)])
        hi = np.log10(angChoice[min(50,minIdx+1)])
        # print lo
        # print hi
        
#     bce = fwd(freq)
#    bce.setspace(nx,ny,dx,dy)
#    bce.setmats(eHS, sHS, nx/2)
#    
#            # pmc = 22.229964825261945
#            # angChoice[minIdx]
#    bce.makeFD(31.2759)
#    bce.point_source(nx/2, ny/2)
#    bce.fwd_solve(0)
    
    # M = bce.nabla2
    
    # matOut.matlab.savemat('bstN', {'M':M})
    
    #do some representative plots
#    plt.figure(1)
#    plt.plot(angChoice, localError)
#    # do some plotting
#    plt.figure(13)
#    plt.subplot(221)
#    plt.imshow((bce.sol[0].real*mask))
#    plt.colorbar()
#    
#    plt.subplot(222)
#    plt.imshow((bbx.real*mask))
#    plt.colorbar()
#    
#    plt.subplot(223)
#    plt.imshow((bce.sol[0].imag*mask))
#    plt.colorbar()
#    
#    plt.subplot(224)
#    plt.imshow((bbx.imag*mask))
#    plt.colorbar()
#    
#    lkl = bce.sol[0]
#    plt.figure(44)
#    plt.plot(np.arange(nx), lkl[49,:].real, np.arange(nx), bbx[49,:].real)
#    
#    plt.figure(4)
#    plt.subplot(121)
#    plt.imshow((bce.sol[0]-bbx).real)
#    plt.colorbar()
#    plt.show()
    
    return angChoice[minIdx]
    

        