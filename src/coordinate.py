'''
Created on Jun 6, 2012
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

Putting together a script to just simply loop through all of these materials and run on the cluster.
Actually, there is no looping. Really, this script is just a wrapper that imports the proper
functions for specifying simulation parameters and calls solveADMM.semiParallel with the
appropriate arguments.

'''


import solveADMM
import sys
import doFolders
import importlib

# trialNs = 100

def main():
    ''' simple main routine '''
    prSpec = importlib.import_module('expts.' + sys.argv[1])

    prSpec.D['outDir'] = doFolders.ensureFolders(prSpec.D, int(sys.argv[2]))
    
    
    gogo = prSpec.getMyVars(int(sys.argv[2]), prSpec.D)
    gogo['ix'] = int(sys.argv[2])    
    solveADMM.semiParallel(**gogo)
    

if __name__ == "__main__":
    main()
