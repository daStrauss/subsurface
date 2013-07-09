'''
Created on Sep 17, 2012
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


D = {'solverType':'sba', 'flavor':'TE', 'numRuns':500, 'expt':'standard', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    bbk = np.int_(np.linspace(10,1600,5))
    romNo, bkg = np.meshgrid(bbk,np.arange(0,100))
    romNo = romNo.flatten()
    bkg = bkg.flatten() 
    
    D['rom'] = romNo[parseNumber]
    D['bkgNo'] = bkg[parseNumber] + 100   
    D['freqs'] = np.array([0.2e3, 6e3, 12e3, 24e3])
    return D
