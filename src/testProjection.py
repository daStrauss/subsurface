'''
Created on Oct 24, 2012

@author: dstrauss
'''


from projection import r5p
import time as tm

def main():
    
    #A = projection.projector()

    xT = 0.5+ 1j*0.7
    yT = 0.2
    zT = -0.3 - 1j*0.2
    xB = 0
    
    strm = tm.time()
    x,y,z = r5p(xT, yT, zT, xB)
    print 'Runtime ' + repr(tm.time()-strm)
    
    print 'x ' + repr(x)
    print 'y ' + repr(y)
    print 'z ' + repr(z)
    
if __name__=="__main__":
    main()