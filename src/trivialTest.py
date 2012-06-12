#!/usr/bin/env python
'''
Created on Jun 11, 2012

@author: dstrauss
'''

import stupidCounter as sc
import time
from mpi4py import MPI

def main():
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    
    ti = time.time()
    sc.busyWork(10000000)
    print repr(rank) + ' of ' + repr(nProc) + ' exec = ' + repr(time.time()-ti)
    
    
    
if __name__ == '__main__':
    main()