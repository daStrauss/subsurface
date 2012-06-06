'''
Created on Jun 6, 2012

@author: dstrauss

Putting together a script to just simply loop through all of these materials and run on the cluster.
'''


import solveADMM

def main():
    ''' simple main routine '''
    
    for ix in range(100):    
        solveADMM.parallel('splitField', rho=1500, xi =2e-4, uBound = 0.05, lmb = 1e-8, bkgNo = (ix+1))


if __name__ == "__main__":
    main()