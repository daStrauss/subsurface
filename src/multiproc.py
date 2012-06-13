from multiprocessing import Pool
import time
import numpy as np
import scipy.sparse
import scipy.sparse.linalg as lin
import time
# from mpi4py import MPI
import scipy.io as spio
# import cProfile

totalIters = 100
N = 2000
# nProc = 2

def zup(inpt):
    # print inpt
    rank = inpt[0]
    nProc = inpt[1]
    print 'rank = ' + repr(rank) + ' nProc = ' + repr(nProc)

    runRange = assgnRange(rank,nProc)
    
    T = np.zeros(totalIters)
    fullT = time.time()
    
    for ix in runRange:
        
        A = scipy.sparse.rand(N,N,0.2)
        A = A + scipy.sparse.spdiags(np.ones(N),0,N,N)
        A = A.tocsc()
        b = np.random.randn(N)
        ti = time.time()
        u = lin.spsolve(A,b)
        T[ix] = time.time()-ti
        
    
    execT = fullT-time.time()
    D = {'times':T,'rng':runRange, 'execT':execT}
    spio.savemat('slvTruePar/tmg'+repr(rank) + '_' + repr(nProc),D)
    print T

def assgnRange(rank, nProc):
    nPer = totalIters/nProc
    leftOver = totalIters - nPer*nProc
    runRange = np.arange(nPer) + rank*nPer
    if leftOver > 0:
        if rank < leftOver:
            runRange = np.append(runRange, (totalIters-1-rank))
    return runRange

def f(r):
    print sum(r)
    
if __name__ == '__main__':
    for nProc in range(1,9):
        p = Pool(nProc)
        a = np.ones((nProc,2),dtype='int')*nProc
        a[:,0] = np.arange(nProc)
        p.map(zup,a)
