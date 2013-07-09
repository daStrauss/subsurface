'''
Created on Jul 18, 2012
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

'''
Created on Jul 16, 2012

@author: dstrauss
'''

import numpy as np


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':613, 'expt':'standard', 'numProcs': 16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    D['bkgNo'] = parseNumber
#    D['rho'] = 0.5e-3
#    D['xi'] = 1e-12
#    D['numProcs'] = 16
#    D['freqs'] = np.array([1e3, 3e3, 13e3, 50e3])  
#    D['inc'] = np.array([75])*np.pi/180
    
    
     
    return D
