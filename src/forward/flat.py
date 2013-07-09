'''
Created on Jun 12, 2012
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

import te
import tm
import te3d

def makeMeA(strg, freq, incAng):
    if strg == 'TE':
        return te.solver(freq, incAng, strg)
    elif strg == 'TM':
        return tm.solver(freq, incAng, strg)
    elif strg == 'TE3D':
        return te3d.solver(freq, incAng, strg)
    else:
        print 'Easy Tiger: ' + repr(strg) + ' aint around'
    
    
