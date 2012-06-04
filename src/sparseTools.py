'''
Created on Jun 3, 2012

@author: dstrauss

I need a few sparse matrix tools

'''

import scipy.sparse as sparse
import numpy as np

def spHcat(A):
    ''' squish matricies given in list A'''
    
    # dims = np.empty(N)
    n = 0
    m = A[0].shape[0]
    # dtx = 0
    I = np.empty(1)
    J = np.empty(1)
    dtx = np.empty(1)
    
    for ix in A:
        I = np.append(I,ix.tocoo().row)
        J = np.append(J,ix.tocoo().col + n)
        dtx = np.append(dtx,ix.tocoo().data)
        n += ix.shape[1]
        if ix.shape[0] != m:
            print 'Major FAIL'
    
    P = sparse.coo_matrix((dtx[1:],(I[1:],J[1:])), shape=(m,n))
    return P
        
def spVcat(A):
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
        
    
        
    
    