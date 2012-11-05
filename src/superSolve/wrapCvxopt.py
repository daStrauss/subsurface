'''
Created on Jul 9, 2012

@author: dstrauss
'''

import cvxopt
from cvxopt import umfpack
import copy
import numpy as np
from cvxopt import lapack
# i guess I don't explicitly check that these are sparse matrices.
# from scipy import sparse

def linsolve(A,b):
    aLocal = A.tocoo()
    AC = cvxopt.spmatrix(aLocal.data.tolist(),aLocal.row.tolist(), aLocal.col.tolist())
    bLocal = cvxopt.matrix(copy.deepcopy(b))
    umfpack.linsolve(AC,bLocal)
    bLocal = np.array(bLocal).flatten()
    return bLocal

def staticSolver(A):
    '''Creates a routine for solving the matrix A --uses UMFPACK underneath'''
    aLocal = A.tocoo()
    AC = cvxopt.spmatrix(aLocal.data.tolist(),aLocal.row.tolist(), aLocal.col.tolist())
    Fs = umfpack.symbolic(AC)
    FA = umfpack.numeric(AC,Fs)
    def Q( b ):
        bLocal = cvxopt.matrix(copy.deepcopy(b))
        umfpack.solve(AC,FA,bLocal)
        bLocal = np.array(bLocal).flatten()
        return bLocal
    
    return Q

def denseSolve(A,b):
    ''' solves an Ax = b matrix system with gesv'''
    if isinstance(A,np.ndarray):
        aLocal = cvxopt.matrix(A)
        bLocal = cvxopt.matrix(b)
        lapack.gesv(aLocal,bLocal)
        return np.array(bLocal).flatten()
    else:
        return linsolve(A,b)

