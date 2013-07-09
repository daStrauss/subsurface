'''
Created on Sep 4, 2012
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


D = {'solverType':'splitField', 'flavor':'TE', 'numRuns':3, 'expt':'spfTest'}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
#    rhos, bkgLocal = np.meshgrid(np.logspace(-4,0,20), range(100))
#    rhos = rhos.flatten()
#    bkgLocal = bkgLocal.flatten()
#    
#    
#    
#    D['bkgSig'] = rhos[parseNumber] 
#    D['bkgNo'] =  bkgLocal[parseNumber] + 100
    D['numProcs'] = 16
    
    if parseNumber == 0:
        D['bkgSig'] = 0.005
        D['bkgNo'] = 100
    elif parseNumber == 1:
        D['bkgSig'] = 0.0048
        D['bkgNo'] = 100
    elif parseNumber == 2:
        # D['bkgSig'] = .5
        D['bkgNo'] = 100
                
    
    return D
