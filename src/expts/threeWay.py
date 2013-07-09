'''
Created on Nov 7, 2012

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
D = {'solverType':'phaseSplit', 'flavor':'TE', 'numRuns':4, 'expt':'testThree'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    if (parseNumber == 0):
        D['freqs'] = np.array([1e3])
        D['numProcs'] = 1
        D['numSensors'] = 2100
    elif (parseNumber == 1):
        D['freqs'] = np.array([1e3, 25e3])
        D['numProcs'] = 2
        D['numSensors'] = 400
    elif (parseNumber == 2):
        D['freqs'] = np.array([25e3])
        D['numProcs'] = 1
        D['numSensors'] = 2100
    elif (parseNumber == 3):
        D['freqs'] = np.linspace(1e3,25e3,6)
        D['numProcs'] = 6
        D['numSensors'] = 400
        
    D['lam'] = 0.0
    D['rho'] = 0.001 
    D['xi'] = 0 
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  100
    D['maxIter'] = 50
    
    return D
