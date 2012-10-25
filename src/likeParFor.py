'''
Created on Oct 24, 2012

@author: dstrauss
'''

from multiprocessing import Pool
import numpy as np
import time

class pp(object):
    def __init__(self,x=0.0+1j*0.0,y=0.0,z=0.0+1j*0.0,xB=0.0+0.0j):
        self.x = x
        self.y = y
        self.z = z
        self.xB = xB
        

def foo(n):
    print 'running ' + repr(n)
    time.sleep(1)
    print 'finishing '
    m = pp()
    m.x = n.x
    m.y = 1.0
    return m
    
    
    
def main():
    pool = Pool(processes=12)
    # output = np.zeros(10)
    
    a = np.random.randn(20)
    b = np.random.randn(20)
    c = np.random.randn(20)
    
    tic = time.time()
    g = [pp(a[itr],b[itr],c[itr]) for itr in range(20)]
    print 'assemble = ' + repr(time.time()-tic)
    
    
    tic = time.time()
    q = pool.map_async(foo, g)
    q.wait()
    
    print repr(q.ready()) + ' ' + repr(time.time()-tic)
    
    # ff = q.get()
    pool.terminate()
    return q,a

if '__name__' == '__main__':
    a = main()