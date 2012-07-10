'''
Created on Jul 9, 2012

@author: dstrauss
'''

import numpy as np
import scipy.sparse
import wrapCvxopt
import scipy.sparse.linalg as lin
import time

def matricles():

    N = 2000
#    A = scipy.sparse.rand(N,N,0.2)
#    A = A + scipy.sparse.spdiags(np.ones(N),0,N,N)
#    A = A.tocoo()
#    AC = cvxopt.spmatrix(A.data.tolist(),A.row.tolist(), A.col.tolist())
#    b = cvxopt.normal(N,1)
    
    for superIx in range(5):
        A = scipy.sparse.rand(N,N,0.2)
        A = A + scipy.sparse.spdiags(np.ones(N),0,N,N)
        A = A.tocsr()
        b = np.random.randn(N)
        
        Q = lin.factorized(A)
        P = wrapCvxopt.staticSolver(A)
        
        for ix in range(1):
            b = np.random.randn(N)
            ti = time.time()
            uOPT = P(b)
            print 'cvx opt time = ' + repr(time.time()-ti)
        
            ti = time.time()
            uSPS = Q(b)
            print 'scipy time  = ' + repr(time.time() - ti)
        
            print np.linalg.norm(uOPT-uSPS)
        
if __name__=='__main__':
    matricles()