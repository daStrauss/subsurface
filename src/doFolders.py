'''
Created on Jun 11, 2012

@author: dstrauss
'''

import sys
import os

def traceBuild(baseTag):
    if not os.path.exists(baseTag):
            os.mkdir(baseTag)
    
    for supIx in range(1):
	#mid = baseTag + '/mat' + repr(supIx) + '/'
	mid = baseTag
        if not os.path.exists(mid):
            os.mkdir(mid)
        
	for ix in range(200):
		outDir = mid + '/trial' + repr(ix) + '/'
		print outDir
                if not os.path.exists(outDir):
                    os.mkdir(outDir)
            
                figDir = outDir + 'Figs/'
                datDir = outDir + 'Data/'
                if not os.path.exists(figDir):
                    os.mkdir(figDir)
                if not os.path.exists(datDir):
                    os.mkdir(datDir)

def main():
    ''' simple main routine '''
    if len(sys.argv) == 1:
        print 'I think you meant to specify one of the following:'
        print 'splitField'
        print 'contrastX'
        print 'sba'
        
    else:
        traceBuild(sys.argv[1] + '/' + sys.argv[2])
     
if __name__ == "__main__":
    main()
