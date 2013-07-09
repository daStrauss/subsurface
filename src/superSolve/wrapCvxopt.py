'''
Created on Jul 9, 2012
Copyright © 2013
The Board of Trustees of The Leland Stanford Junior University.
All Rights Reserved

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

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

def createSymbolic(A):
    ''' returns a symbolic factorization object for later reuse'''
    
    s = A.shape
    aLocal = A.tocoo()
    AC = cvxopt.spmatrix(aLocal.data.tolist(),aLocal.row.tolist(), aLocal.col.tolist(),s)
    Fs = umfpack.symbolic(AC)
    return Fs

def solveNumeric(A,b, Fs):
    ''' given a static Fs, or symbolic factorization of the matrix A, performs the numeric part '''
    aLocal = A.tocoo()
    s = A.shape

    AC = cvxopt.spmatrix(aLocal.data.tolist(),aLocal.row.tolist(), aLocal.col.tolist(),s)
    # Fs = umfpack.symbolic(AC)
    FA = umfpack.numeric(AC,Fs)
    bLocal = cvxopt.matrix(copy.deepcopy(b))
    umfpack.solve(AC,FA,bLocal)
    bLocal = np.array(bLocal).flatten()
    return bLocal
    

def denseSolve(A,b):
    ''' solves an Ax = b matrix system with gesv'''
    if isinstance(A,np.ndarray):
        aLocal = cvxopt.matrix(A)
        bLocal = cvxopt.matrix(b)
        lapack.gesv(aLocal,bLocal)
        return np.array(bLocal).flatten()
    else:
        return linsolve(A,b)

