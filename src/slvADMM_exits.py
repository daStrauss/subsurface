'''
Created on Jun 15, 2012

@author: dstrauss
'''

def serial(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1):

    # rho = 1500.0
    # ho = 1e-3
    # xi = 2e-3
    # lmb = 0.0
    
    print 'xi = ' + repr(xi) + ' rho = ' + repr(rho)

    freq = np.array([1e3, 1e4])
    incAng = np.array([45.0, 45.0])*np.pi/180.0
    ranks = np.arange(np.size(freq))
    
    S = delegator(solverType, freq, incAng, ranks)
    # S = smallProj(S)
    S,pTrue = bigProj(S,'basic',bkgNo)
    N = np.size(S)
    
    for F in S:
        uHat = F.Ms*(F.sol[1].flatten())
        ti = time.time()
        F.initOpt(uHat,rho,xi,uBound, lmb)
        print 'initalization time ' + repr(time.time()-ti)
    
    # P = np.zeros(80*25)
    P = np.zeros(S[0].nRx*S[0].nRy)
    
    resid = np.zeros(100)
    
    # io.savemat('bkg', {'u':S[0].sol[1]})
    # io.savemat('sigm', {'sigMap':S[0].sigmap[1]})
    
    for itNo in range(100):
        
        ti = time.time()
        for ix in range(N):
            # run each individual update
            S[ix].runOpt(P)
        
        print 'iter no ' + repr(itNo) + ' exec time ' + repr(time.time()-ti)
        # aggregate over all
        if solverType == 'sba':
            P += S[0].aggregatorSerial(S)
        else:
            P = S[0].aggregatorSerial(S)
            
        resid[itNo] = np.linalg.norm(P-pTrue)

    S[0].plotSerial(S, P, resid)


def parallel(solverType, rho=1e-3, xi=2e-3, uBound=0.05, lmb=0, bkgNo=1, outDir='basic'):
    '''Parallel solver -- i.e. has MPI calls'''
    
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    nProc = comm.Get_size()
    timeFull = time.time()
    
    fout = open('notes'+repr(rank) + '_' +repr(bkgNo) + '.nts', 'w')
    
    fout.write('xi ' + repr(xi) + ' rho = ' + repr(rho) + '\n')
        
    #  
    allFreq = np.array([1e3, 3e3, 13e3, 50e3])
    allIncAng = np.ones(allFreq.shape)*45*np.pi/180.0
    # allRanks = np.arange(np.size(freq))
    
    S = delegator(solverType, [allFreq[rank]], [allIncAng[rank]], [rank])
    S,pTrue = bigProj(S, outDir, bkgNo)
    
    # de reference so that I don't continuously have to work with lists in parallel mode
    S = S[0]
    
    uHat = S.Ms*(S.sol[1].flatten())
    
    ti = time.time() 
    S.initOpt(uHat, rho, xi, uBound, lmb)
    fout.write('initialization time ' + repr(time.time()-ti) + '\n')
    
    P = np.zeros(S.nRx*S.nRy)
    resid = np.zeros(50)
    
    for itNo in range(50):
        ti = time.time()        
        S.runOpt(P)
        
        # i don't think i can get around this!
        if solverType == 'sba':
            P += S.aggregatorParallel(comm)
        else:
            P = S.aggregatorParallel(comm)
            
        resid[itNo] = np.linalg.norm(P-pTrue)
        fout.write('iter no ' + repr(itNo) + ' exec time = ' + repr(time.time()-ti) + ' rank ' + repr(comm.Get_rank()) +'\n')
        
    # do some plotting        
    S.plotParallel(P,resid,rank)
    S.writeOut()
    
    if rank == 0:
        D = {'Pfinal':P.reshape(S.nRx,S.nRy), 'nProc':nProc, 'resid':resid}
        spio.savemat('pout', D)
        
    fout.write('Solve time = ' + repr(time.time()-timeFull) + '\n')
    fout.close()