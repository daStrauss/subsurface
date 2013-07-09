'''
Created on Oct 5, 2012

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
import scipy.io as spio

# F = spio.loadmat('incCondGo.mat')
# numRuns = F['goTo'].shape[0]


D = {'solverType':'contrastX', 'flavor':'TE', 'numRuns':4800, 'expt':'noSense', 'numProcs':16}


def getMyVars(parseNumber, D):
    '''routine to return the parameters to test at the current iteration.'''
    # noFreqs,noPhis,bkg = np.meshgrid(range(1,7), range(1,7), range(100))
    noFreqs,noPhis,bkg = np.mgrid[1:7,0:16,0:50]
    noFreqs = noFreqs.flatten()
    noPhis = noPhis.flatten() 
    bkg = bkg.flatten()
    
    allfrq = np.arange(1,80,5)
        
    
    D['freqs'] = np.round(np.logspace(np.log10(1000), np.log10(50000), noFreqs[parseNumber]))
    # D['inc'] = [75.0*np.pi/180.0] -- use defaults!
    D['numSensors'] = allfrq[noPhis[parseNumber]]
    
    D['bkgNo'] = bkg[parseNumber]+100;
    D['numProcs'] = len(D['freqs'])*4
            
        
    return D
