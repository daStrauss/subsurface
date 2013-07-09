'''
Created on Jan 2, 2013

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

D = {'solverType':'sba', 'flavor':'TE', 'numRuns':2000, 'expt':'bringDa', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    
    snr,bkg = np.meshgrid(np.logspace(-5,0,20),range(1000))
        
    snr = snr.flatten()
    bkg = bkg.flatten()                 

    D['relNoi'] = snr[parseNumber]        
    D['bkgNo'] = bkg[parseNumber]+100

    
        
    return D
