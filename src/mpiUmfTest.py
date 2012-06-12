'''
Created on Jun 7, 2012

@author: dstrauss
'''

import numpy
import scipy.sparse
import scipy.sparse.linalg as lin
import time
from mpi4py import MPI
import cProfile

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    # nProc = comm.Get_size()
    # rank = 5
    A = numpy.random.randn(5000,5000)
    A[A<0.5] = 0
    
    M = scipy.sparse.csc_matrix(A)
    
    f = open('./notes' + repr(rank) + '.txt', 'w')
    
    
    for ix in range(20):
        ti = time.time()
        Z = lin.factorized(M)
        f.write('rank ' + repr(rank) + ' iter ' + repr(ix) + ' exec time = ' + repr(time.time()-ti) + '\n')
    
    f.close()
        
if __name__ == "__main__":
    cProfile.run('main()')