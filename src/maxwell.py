'''
Created on May 23, 2012

@author: dstrauss
'''
import numpy as np
import scipy as sp
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
        
        
        