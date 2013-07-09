'''
Created on Jan 9, 2013
Copyright © 2013
The Board of Trustees of The Leland Stanford Junior University.
All Rights Reserved

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at
       http://www.apache.org/licenses/LICENSE-2.0
Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

@author: dstrauss
'''

import numpy as np
D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':100, 'expt':'goBig', 'numProcs':16}



def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    

        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), 100))
    D['inc'] = [75.0*np.pi/180.0] # one freq.
    D['numSensors'] = 20
    D['bkgNo'] = parseNumber+100;
    D['numProcs'] = 100
    D['rho'] = 1e-3
    D['xi'] = 1e-12

    
            
        
    return D
