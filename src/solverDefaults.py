'''
Created on Jun 30, 2012

@author: dstrauss
'''

def getDefaults(solver,flavor):
    ''' just a stub function that will allow me to produce a dictionary that 
    gives default values'''
    
    if (solver == 'splitField') & (flavor == 'TE'):
        D = {'rho':148, 'xi':1e-3, 'lmb':1e-8, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/'}
        return D
    
    elif (solver =='splitField') & (flavor == 'TM'):
        D = {'rho':1 xi':1e-3, 'lmb':1e-8, 'uBound':0.05, 'bkgNo':1, 'outDir':'testing/'}