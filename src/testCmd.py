'''
Created on Jun 8, 2012

@author: dstrauss
'''
import sys

def zipps(frg='def'):
    print 'you wrote: ' + frg
    print 'boog ' + \
        'does split line work?'
    
if __name__ == '__main__':
    if len(sys.argv)   > 1:
        
        zipps(sys.argv[1])
    else:
        zipps()
        