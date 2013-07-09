'''
Created on Nov 28, 2012

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


D = {'solverType':'phaseSplitSingle', 'flavor':'TE', 'numRuns':960, 'expt':'paramPSS', 'numProcs':1}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    rhos, xis, fr, bkgLocal = np.meshgrid(np.logspace(-4,4,8), np.linspace(0,1,8), np.array([1e3, 1e4, 1e5]), range(5))
    rhos = rhos.flatten()
    xis = xis.flatten()
    bkgLocal = bkgLocal.flatten()
    fr = fr.flatten()
        
    D['freqs'] = np.array(fr[parseNumber])
    D['numProcs'] = 1
    D['numSensors'] = 3100

    D['lam'] = 0
#    D['rho'] = 1e-3    
    D['inc'] = np.array([75*np.pi/180])
    D['maxIter'] = 200
    
    D['rho'] = rhos[parseNumber] 
    D['xi'] = xis[parseNumber]
    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    return D
