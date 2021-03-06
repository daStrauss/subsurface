'''
Created on Jun 3, 2012
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

I need a few sparse matrix tools

'''

import scipy.sparse as sparse
import numpy as np

def hCat(A):
    ''' squish matricies given in list A'''
    
    # dims = np.empty(N)
    n = 0
    m = A[0].shape[0]
    # dtx = 0
    I = np.empty(1)
    J = np.empty(1)
    dtx = np.empty(1)
    
    for ix in A:
        if not sparse.isspmatrix(ix):
            ix = sparse.coo_matrix(ix)
        
        I = np.append(I,ix.tocoo().row)
        J = np.append(J,ix.tocoo().col + n)
        dtx = np.append(dtx,ix.tocoo().data)
        n += ix.shape[1]
        if ix.shape[0] != m:
            print 'Major FAIL'
    
    P = sparse.coo_matrix((dtx[1:],(I[1:],J[1:])), shape=(m,n))
    return P
        
def vCat(A):
    '''Stack matrices given in list A '''
    n = A[0].shape[1]
    m = 0
    I = np.empty(1)
    J = np.empty(1)
    dtx = np.empty(1)
    
    for ix in A:
        lcl = ix.tocoo()
        I = np.append(I,lcl.row + m)
        J = np.append(J,lcl.col)
        dtx = np.append(dtx,lcl.data)
        m += ix.shape[0]
        if ix.shape[1] != n:
            print 'Major FAIL'
        
    P = sparse.coo_matrix((dtx[1:],(I[1:],J[1:])), shape=(m,n))
    return P 
        
    
def smartX(A,B):
    if sparse.issparse(A) | sparse.issparse(B):
        return A*B
    elif isinstance(A,np.ndarray) & isinstance(B,np.ndarray):
        return np.dot(A,B)
    else:
        print 'Take your chances'
        return A*B
        
    
    
