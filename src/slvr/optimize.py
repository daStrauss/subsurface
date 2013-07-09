'''
Created on Jun 14, 2012

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

abstract definition for specifying any given optimization object.
Not a formal requirement, since Python can be pretty loose with
abstract objects. 
'''

import forward.flat

class optimizer(object):
    '''base class for optimization routines'''
    
    def __init__(self,freq,incAng,flavor):
        self.fwd = forward.flat.makeMeA(flavor, freq, incAng)
    
    
