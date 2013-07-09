'''
Created on Oct 17, 2012
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
D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':2, 'expt':'testSolver'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
#    D['freqs'] = np.array([1e3])
#    D['inc'] = np.array([75*np.pi/180])
#    D['bkgNo'] =  100
#    D['numProcs'] = 1
#    D['maxIter'] = 100
#    D['numSensors'] = 3010
#    D['rho'] = 0.10
    if parseNumber == 0:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([75*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 1e-3
        D['xi'] = 1e-12
    elif parseNumber == 1:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([75*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 5e-4
        D['xi'] = 1e-12
        
    elif parseNumber == 2:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([45*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 1e-3
        D['xi'] = 1e-12
    elif parseNumber == 3:
        D['freqs'] = np.array([1000, 3684, 13572, 50000])
        D['inc'] = np.array([45*np.pi/180])
        D['bkgNo'] =  100
        D['numProcs'] = 4
        D['maxIter'] = 200
        D['numSensors'] = 61
        D['rho'] = 5e-4
        D['xi'] = 1e-12
    
    return D
