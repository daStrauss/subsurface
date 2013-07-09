'''
Created on Nov 8, 2012

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


D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':245, 'expt':'newParams2', 'numProcs':6}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis, bkgLocal = np.meshgrid(np.logspace(-6,-3,7), np.logspace(-14,-11,7), range(5))
    rhos = rhos.flatten()
    xis = xis.flatten()
    bkgLocal = bkgLocal.flatten()
    
        
    D['freqs'] = np.linspace(1e3,25e3,6)
    D['numProcs'] = 6
    D['numSensors'] = 400

    D['lam'] = 0
#    D['rho'] = 1e-3    
    D['inc'] = np.array([75*np.pi/180])
#    D['maxIter'] = 200

    
    
    D['rho'] = rhos[parseNumber] 
    D['xi'] = xis[parseNumber]
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    return D
