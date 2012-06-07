'''
Created on Jun 7, 2012

@author: dstrauss
'''

import numpy
import scipy.sparse
import scipy.sparse.linalg as lin
import time


def main():
    A = numpy.random.randn(1000,1000)
    A[A<0] = 0
    
    M = scipy.sparse.csc_matrix(A)
    
    for ix in range(10000):
        ti = time.time()
        
        Z = lin.factorized(M)
        print repr(ix) + ' exec time = ' + repr(time.time()-ti)
        
if __name__ == "__main__":
    main()