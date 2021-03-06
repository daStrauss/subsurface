'''
Created on Nov 14, 2012

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

Simple routine to do optimization over 3D! 

'''


import numpy as np
D = {'solverType':'splitField', 'flavor':'TE3D', 'numRuns':4, 'expt':'huge', 'numProcs':1}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    
    
    # D['rho'] = 0.00001
    # D['xi'] = 1e-9
    D['freqs'] = np.array([1e3])
    D['inc'] = np.array([75*np.pi/180])
    D['bkgNo'] =  0
    D['numProcs'] = 1
    D['reg'] = 1e-6
    if parseNumber == 0:
        D['maxIter'] = 1
    elif parseNumber == 1:
        D['maxIter'] = 200
    elif parseNumber == 2:
        D['maxIter'] = 10
    elif parseNumber == 3:
        D['maxIter'] = 10
        D['reg'] = 0.0   
        
    elif parseNumber == 4:
        D['bkgSig'] = 1e-5
        D['maxIter'] = 30
        D['freqs'] = np.array([1e5])
        D['reg'] = 0.0 
        D['rho'] = 100.0
        
        
    D['lmb'] = 0e-5
    
    return D
