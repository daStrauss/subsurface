'''
Created on Jun 21, 2012

@author: dstrauss
'''

from mpi4py import MPI
import time
N = 1e8


def main():
    '''simple function to wait and do nothing '''
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    name = MPI.Get_processor_name()
    
    
    print 'I am ' + name + ' number ' + repr(rank) + ' of ' + repr(nProc)
    n = N
    ti = time.time()
    while n >= 0:
        n -= 1
    print time.time() - ti
        
        
if __name__=='__main__':
    main()