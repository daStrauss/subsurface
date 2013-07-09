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
        
        
        
        
        for ix in range(1):
            b = np.random.randn(N)
            ti = time.time()
            P = wrapCvxopt.staticSolver(A)
            uOPT = P(b)
            print 'cvx opt time = ' + repr(time.time()-ti)
        
            ti = time.time()
            Q = lin.factorized(A)
            uSPS = Q(b)
            print 'scipy time  = ' + repr(time.time() - ti)
        
            print np.linalg.norm(uOPT-uSPS)

def reNumber():
    N = 2000
    A = scipy.sparse.rand(N,N,0.2) + scipy.sparse.eye(N,N)
    b = np.random.randn(N)
    
    aLocal = A.tocoo()
    
    I = aLocal.row.tolist()
    J = aLocal.col.tolist()
    D = aLocal.data.tolist()
    
    print len(D)
    # print A.shape
    alltm = time.time()
    Fs = wrapCvxopt.createSymbolic(A)
    solA = wrapCvxopt.solveNumeric(A, b, Fs)
    # return A,Fs
    # print A.shape
    print 'symb and num ' + repr(time.time()-alltm)
    
    lnf = time.time()
    Q = lin.factorized(A)
    uSPS = Q(b)
    print 'umfpack + scipy ' + repr(time.time()-lnf)
    

    M = scipy.sparse.coo_matrix((np.random.randn(len(D)),(I,J)))
    numTim = time.time()
    solB = wrapCvxopt.solveNumeric(M, b, Fs)
    print Fs
    print 'num only ' + repr(time.time()-numTim)
    
    lnf = time.time()
    Q = lin.factorized(M)
    uSPS = Q(b)
    print 'umfpack + scipy ' + repr(time.time()-lnf)
    
        
if __name__=='__main__':
    matricles()
