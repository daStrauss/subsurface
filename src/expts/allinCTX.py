'''
Created on Jan 3, 2013

@author: dstrauss
'''

import numpy as np
D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':100, 'expt':'noSense', 'numProcs':16}



def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), 100))
    D['inc'] = [75.0*np.pi/180.0] # one freq.
    D['numSensors'] = 20
    D['bkgNo'] = parseNumber+100;
    D['numProcs'] = 100
            
        
    return D