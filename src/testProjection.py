'''
Created on Oct 24, 2012

@author: dstrauss
'''


import projection
import time as tm

def main():
    
    A = projection.projector()

    xT = 0.5+ 1j*0.7
    yT = 0.2
    zT = -0.3 - 1j*0.2
    xB = 0
    
    str = tm.time()
    x,y,z = A.r5p(xT, yT, zT, xB)
    print 'Runtime ' + repr(tm.time()-str)
    
    print 'x ' + repr(x)
    print 'y ' + repr(y)
    print 'z ' + repr(z)
    
if __name__=="__main__":
    main()