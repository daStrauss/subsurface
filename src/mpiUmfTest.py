'''
Created on Jun 7, 2012

@author: dstrauss
'''

import numpy as np
import scipy.sparse
import scipy.sparse.linalg as lin
import time
from mpi4py import MPI
import scipy.io as spio
# import cProfile

totalIters = 20

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    # rank = 5
    
    runRange = assgnRange(rank,nProc)
    
    T = np.zeros(totalIters)
    fullT = time.time()
    
    for ix in runRange:
        ti = time.time()
        A = scipy.sparse.rand(3000,3000,0.2).tocsc()
        b = np.random.randn(3000)
        u = lin.spsolve(A,b)
        T[ix] = time.time()-ti
        
    
    execT = fullT-time.time()
    D = {'times':T,'rng':runRange, 'execT':execT}
    spio.savemat('tmg/tmg'+repr(rank) + '_' + repr(nProc),D)
    

def assgnRange(rank, nProc):
    nPer = totalIters/nProc
    leftOver = totalIters - nPer*nProc
    runRange = np.arange(nPer) + rank*nPer
    if leftOver > 0:
        if rank < leftOver:
            runRange = np.append(runRange, (totalIters-1-rank))
    return runRange
    
    
    
if __name__ == "__main__":
    main()