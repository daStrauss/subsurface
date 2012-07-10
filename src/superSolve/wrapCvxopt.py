'''
Created on Jul 9, 2012

@author: dstrauss
'''

import cvxopt
from cvxopt import umfpack
import copy
import numpy as np


def linsolve(A,b):
    aLocal = A.tocoo()
    AC = cvxopt.spmatrix(aLocal.data.tolist(),aLocal.row.tolist(), aLocal.col.tolist())
    bLocal = cvxopt.matrix(copy.deepcopy(b))
    umfpack.linsolve(AC,bLocal)
    bLocal = np.array(bLocal).flatten()
    return bLocal

def staticSolver(A):
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

