'''
Created on Jun 30, 2012

@author: dstrauss
'''
import numpy as np

def getDefaults(solverType, flavor, customDefs):
    ''' just a stub function that will allow me to produce a dictionary that 
    gives default values, also populates with new preferences.'''
    
    if (solverType == 'splitField') & (flavor == 'TE'):
        D = {'rho':2320.27, 'xi':9.811e-4,'lmb':0.0}
   
    elif (solverType =='splitField') & (flavor == 'TM'):
        D = {'rho':0.019307, 'xi':1.3895e-3, 'lmb':0.0}
    
    elif (solverType =='sba'):
        D = {'rho':0.005, 'xi':0.9, 'lmb':0} 
    
    elif (solverType =='contrastX') & (flavor =='TE'):
        D = {'rho':1e-3, 'xi':2e-3, 'lmb':0}
    
    elif (solverType =='biconvex') & (flavor=='TE'): 
        D = {'rho':2069.1, 'xi':1e-5, 'lmb':0}
        
    elif (solverType =='contrastX') & (flavor =='TM'):
        D = {'rho':1e-3, 'xi':2e-3, 'lmb':0}
        
    elif (solverType =='splitField') & (flavor =='both'):
        D = {'rho':1e-3, 'xi':1e-3,'lmb':0}
        
    elif (solverType == 'contrastSoftX') & (flavor=='TE'):
        D = {'rho':1e-3, 'xi':1e-3,'lmb':0}
        
    elif (solverType == 'projection') & (flavor=='TE'):
        D = {'rho':1e-3, 'xi':1e-3, 'lmb':0}
    elif (solverType == 'middleMan') & (flavor=='TE'):
        D = {'rho':1e-3, 'xi':1e-3, 'lmb':0}
    
    else:
        print 'somehow you did not specify a valid combination'
        

    
    D['bkgNo'] = 1
    D['outDir'] = 'testing/'
    D['uBound'] = 0.05
    D['rom'] = False
    D['freqs'] = np.array([1e3, 3e3, 13e3, 50e3])  
    D['inc'] = np.array([75, -75, 45, -45])*np.pi/180
    D['numProcs'] = 16
    D['maxIter'] = 1000
    D['bkgSig'] = 0.005;
    D['flavor'] = flavor
    D['numSensors'] = 10
    # D['flavor'] = 'TE' # do this by default, although it should be set explicitly
    
    for key,value in customDefs.iteritems():
        D[key] = value
      
    return D
        