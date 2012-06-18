'''
Created on Jun 18, 2012

@author: dstrauss
'''

import xml.etree.ElementTree as xml
import subprocess
import time
 
    
def waitForExit(jobName):
    ''' Routine to check and see if a particular job is still running, if not, return '''
    doExit = False
    while doExit == False:
        xmlOut = subprocess.check_output('qstat ' + jobName + ' -x')
        P = xml.parse(xmlOut)
        for z in P.iter():
            if z.tag == 'job_state':
                if z.text == 'R':
                    print 'Still Running'
                elif z.text == 'C':
                    doExit = True
                    print 'Finished ' + jobName
        
        time.sleep(5)
        
    return 1


if __name__=='__main__':
    jobName = '33591.nansen'
    waitForExit(jobName)