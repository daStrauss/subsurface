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
        
    def foo(self):
        ''' should do something simple '''
        return self.x*self.y
        

def foo(m):
    print 'running ' + repr(m[0])
    time.sleep(1)
    print 'finishing '
    return m[0]*m[1],m[0]+m[1]
    
    
    
def main():
    pool = Pool(processes=12)
    # output = np.zeros(10)
    
    a = np.arange(20)
    b = np.arange(20) + 100
    c = np.random.randn(20)
    
    tic = time.time()
    g = [pp(a[itr],b[itr],c[itr]) for itr in range(20)]
    print 'assemble = ' + repr(time.time()-tic)
    
    
    tic = time.time()
    q = pool.map_async(foo, zip(a,b))
    q.wait()
    
    print repr(q.ready()) + ' ' + repr(time.time()-tic)
    t = q.get()
    
    ao = np.zeros(20)
    bo = np.zeros(20)
    
    for ixr,res in enumerate(t):
        ao[ixr] = res[0]
        bo[ixr] = res[1]
        
    # ff = q.get()
    pool.terminate()
    return ao,bo

if '__name__' == '__main__':
    a = main()