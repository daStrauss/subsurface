'''
Created on Jun 30, 2012

@author: dstrauss
'''
import numpy as np

def getDefaults(solver,flavor,**kwargs):
    ''' just a stub function that will allow me to produce a dictionary that 
    gives default values, also populates with new preferences.'''
    
    if (solver == 'splitField') & (flavor == 'TE'):
        D = {'rho':148, 'xi':1e-3,\
            'lmb':1e-8, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/', \
            'rom':False}

    
    elif (solver =='splitField') & (flavor == 'TM'):
        D = {'rho':0.019307, 'xi':1.3895e-3, \
             'lmb':1e-8, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/',\
             'rom':False}
    
    elif (solver=='sba'):
        D = {'rho':0.005, 'xi':0.9, \
             'lmb':0, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/',\
             'rom':False} 
    
    elif (solver=='contrastX') & (flavor =='TE'):
        D = {'rho':1e-3, 'xi':2e-3, \
             'lmb':0, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/',\
             'rom':False}
    
    elif (solver=='biconvex') & (flavor=='TE'): 
        D = {'rho':0.001, 'xi':1e-5, \
             'lmb':0, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/', \
             'rom':False}
    
    else:
        print 'somehow you did not specify a valid combination'
        
    for key,value in kwargs.iteritems():
        D[key] = value
    
    D['freqs'] = np.array([1e3, 3e3, 13e3, 50e3])  
    D['inc'] = np.array([75, -75, 45, -45])*np.pi/180
      
    return D
        