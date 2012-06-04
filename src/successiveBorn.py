'''
Created on Jun 3, 2012

@author: dstrauss

Implementation of the successive born method

'''

import numpy as np
from maxwell import twoDim
import scipy.sparse as sparse
import scipy.sparse.linalg as lin
from mpi4py import MPI

class sla(twoDim):
    ''' class for implementing successive linear approximation for the
    subsurface imaging problem '''
    def initOpt(self, rho, alpha, uHat):
        self.rho = rho
        self.alpha = alpha
        self.uHat = uHat
        

    
    def runOpt(self, P):
        
        