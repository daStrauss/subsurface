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
import superSolve.wrapCvxopt
import time
from scipy import linalg

class pmlList(object):
    freq = np.ndarray(0)
    dx = np.ndarray(0)
    pmc = np.ndarray(0)

class fwd(object):
    epso = 8.854e-12
    muo = 4*np.pi*1e-7
    c = np.sqrt(1/(muo*epso))
    
    def __init__(self, freq, incAng, flavor, rank=0, azAng=0.0):   
        '''initialize object ''' 
        self.f = freq
        self.w = 2*np.pi*freq
        self.l = self.c/self.f
        self.incAng = incAng
        self.azAng = azAng # add a spot for the azimuth angle in 3D
        self.gogo = ['', '']
        self.rom = False
        self.flavor = flavor
        
    def initBig(self, p, bkgSig=0.005, numSensors=10):
        '''Create a "big" (nx=ny=199) style problem with some basic background parameters '''
        self.setspace(199,199,5.0,5.0)
        self.setmats(1,bkgSig,199/2)
        self.setOperators()
        self.makeGradOper()
        self.setMs(numSensors)
        self.setMd([60,140], [70,95])
        self.sigmap[1] += self.Md.T*p
        
        self.planeWave()
        self.fwd_solve(0)
        self.fwd_solve(1)
        
    def initSmall(self, p, bkg=0.005, numSensors=10):
        '''Create a "big" (nx=ny=199) style problem with some basic background parameters '''
        self.setspace(99,99,5.0,5.0)
        self.setmats(1,bkg,99/2)
        self.setOperators()
        self.makeGradOper()
        self.setMs(numSensors)
        self.setMd([30, 70], [35, 45])
        self.sigmap[1] += self.Md.T*np.ones(40*10)*0.01 # don't worry about specific shapes.        
        self.planeWave()
        self.fwd_solve(0)
        self.fwd_solve(1)
        
    def init3Dfeas(self, p, bkg=0.005, numSensors=10):
        '''Create a "big" (nx=ny=199) style problem with some basic background parameters '''
        self.setspace(31,31,31,5.0,5.0,5.0)
        self.setmats(1,bkg,31/2)
        self.setOperators()
        self.makeGradOper()
        self.setMs(numSensors)
        ''' for 31,31,31 -- use '''
        self.setMd([9,23], [8,15], [9,23])
        # self.sigmap[1] += self.Md.T*np.ones(1372)*0.01 # don't worry about specific shapes. 
        ''' for 21, 21,21 -- use '''
        #self.setMd([7,14], [6,9],[7,14])
        
        ''' for 41,41,41 '''
#        self.setMd([12,30], [12,19], [12,30])
        
        
        self.sigmap[1] += self.Md.T*p # don't worry about specific shapes. 
        self.planeWave()
        self.fwd_solve(0)
        self.fwd_solve(1)
         
    def setspace(self, nx,ny,dx,dy):
        self.nx = nx # number x points
        self.ny = ny # number y points
        self.N = nx*ny # super size of space
        self.dx = dx # delta x
        self.dy = dy # delta y

        self.npml = min(10,round((nx+2.0)/10)) # size of pml borders
    
    def setmats(self, eHSr, sHS, div):
        """Routine to setup and allocate the materials - epsilon and sigma. 
        Only sigma is implemented for the TM mode. Settings taken for ON Grid Locations"""
        pass
        
    def getk(self,ind):
        '''Routine to assemble diagonal matrix that adds in the materials for the entire space.
        Indexed by ind for the two spaces held by the model'''
        pass    

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
            lkt.dx = np.append(lkt.dx, self.dx)
            pickle.dump(lkt, open('pmlLib.p', 'wb'))
            # print pmc
            return pmc
            
    def setOperators(self):
        ''' An wraper for makeFinteDiferences that will select the right PML parameters
        '''
        # self.makeFD(1886.7678)
        self.makeFD(self.getPMLparm())    
            
        
    def makeFD(self, pmc):
        ''''A routine that will define the following operators:
            d1: derivative mapping from on grid to half grid
            d2: derivative mapping from half grid to on grid
            po: pml scaling for on grid
            ph: pml scaling for half grid
        '''
            
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
                 
    def pointSource(self, x,y):
        ''' This routine puts a point source in the z direction, whichever field that may be'''
        print 'Not yet implemented pointSource'
        pass    

    def fwd_solve(self, ind):
        '''Does the clean solve for the given index. The factorization is not cached'''      
        strm = time.time()
        self.gogo[ind] = lin.factorized(sparse.csc_matrix(self.nabla2+self.getk(ind)) )
        print 'factor time = ' + repr(time.time()-strm)
        # self.gogo[ind] = superSolve.wrapCvxopt.staticSolver(self.nabla2+self.getk(ind))
        strm = time.time()
        self.sol[ind] = self.gogo[ind](self.rhs.flatten())
        print 'sol time ' + repr(ind) + ' time = ' + repr(time.time()-strm)
        
    def parseFields(self,u):
        '''method to reshape and return according to internal dimensions '''
        print 'not yet implemented parseFields'
        pass

    def setMs(self, nSensors=10):
        '''Tell me the number of sensors, and I will distribute them equally across the surface'''
        pass
         
    def setMd(self, xrng, yrng):
        '''Tell me the xrange and the yrange and I'll make selector'''  
        pass
        
    def plotSol(self,ind):
        pass

    def writeOut(self,ind):
        '''write out the juicy parts in matlab format'''
        D = {'f':self.f, 'angle':self.incAng, 'sigMat':self.simap[ind], 'fld':self.sol[ind]}
        spio.savemat('maxPrint' + repr(self.rank), D)
        
    def buildROM(self,nBases, force=True):
        ''' Extending to be able to do the reduced order model:
        buildRom(nBases, force) : nBases is the number of modes to keep, force means should I really do it
        or look for something in the files.'''
        self.rom = nBases
        ti = time.time()
        if nBases == self.N:
            self.Phi = sparse.eye(self.N,self.N,dtype='complex128')
            print 'In which we cheat. time = ' + repr(time.time()-ti)
            return
            
        if not force:
            try:
                
                P = pickle.load(open('teCache.p','rb'))
                self.Phi = P
                print 'loaded ROM from cache'
            except:
                print 'cannot find teCache'
                force = True
                 
        
        if force:
            # do all of the points in the Md region
            M = self.nRx*self.nRy
            X = np.zeros((self.N,M+1),dtype='complex128')
            X[:,0] = self.sol[0]
            
            for ix in range(M):
                print 'computing derivative solution ' + repr(ix) + ' of ' + repr(M)
                p = np.zeros(M,dtype='complex128')
                p[ix] = 1j/self.w
                X[:,ix+1] = self.gogo[0](self.Md.T*p)
            
            print 'Computing SVD'
            u,s,v = linalg.svd(X,full_matrices=False)
            
            self.Phi = u[:,:nBases]
            # print self.Phi.dtype
            # pickle.dump(self.Phi, open('teCache.p','wb')) 
            
            
        print 'ROM build time = ' + repr(time.time()-ti)
                     
    
def findBestAng(freq,dx):
    '''for a particular frequency, find the best PML complex angle 
    Universal for both te, tm models, but calculates based on the TE model
    '''
    import te 
    print 'blast -- recalculating proper PML for f = ' + repr(freq)
    # set some internal parameters
    nx = 99
    ny = nx
    dy = dx
    
    # this means free space
    eHS = 1.0
    sHS = 0
    
    # need basic constants
    epso = 8.854e-12
    muo = 4.0*np.pi*1e-7
    c = 1/np.sqrt(muo*epso)
    
    # to calculate the Greens function, the true spatial locations are needed
    x = np.array(range(nx))*dx - dx*(nx/2)
    Y,X = np.meshgrid(x, x);
    dist = np.sqrt(Y**2+X**2) 
    
    k = (2*np.pi)/(c/freq)
    bbx = dx*dx*(1j/4)*spec.hankel1(0,k*dist)
    mdpt = nx/2
    # the center contains a singularity, but we'll  be forgiving and finite
    bbx[mdpt,mdpt] = 0.25*bbx[mdpt+1,mdpt] + 0.25*bbx[mdpt-1,mdpt] + 0.25*bbx[mdpt,mdpt+1] + 0.25*bbx[mdpt,mdpt-1]
    
    # create a mask so that PML effects don't adversely affect the quality of the solutions
    mask = np.ones(bbx.shape)
    mask[:15,:] =0
    mask[(nx-15):(nx),:] = 0
    mask[:,:15]=0
    mask[:,(nx-15):(nx)] = 0
    
    lo = 0
    hi = 5
    
    # do two loops over different ranges to get a more precise estimate
    for tune in range(2):
        localError = np.zeros(50)
        angChoice = np.logspace(lo,hi,50)
    
        for i in range(50):    
            # Create a single new object on frequency
            # use the te routine for this
            bce = te.solver(freq,0.0,'TE')
            bce.setspace(nx,ny,dx,dy)
            bce.setmats(eHS, sHS, nx/2)
            
            bce.makeFD(angChoice[i])
            bce.makeGradOper()
            bce.pointSource(49, 49)
            bce.fwd_solve(0)
            u = bce.parseFields(bce.sol[0])
            
            localError[i] = np.linalg.norm((u[0] - bbx)*mask,'fro')
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
    

        