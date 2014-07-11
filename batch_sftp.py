__author__ = 'brad'

''' Mass ftp script for vtp streamline files.'''

#http://docs.paramiko.org/en/1.13/api/sftp.html

__author__ = 'brad'

from ftplib import FTP

import os
import sys
import time
from subprocess import *
import paramiko

def ftpConnect():
    
    global PASSWORD
    
    host ='sundance.soe.ucsc.edu'
    port = 22

    username = 'behollis'

    transport = paramiko.Transport((host, port))
    transport.connect(username = username, password = PASSWORD)

    sftp = paramiko.SFTPClient.from_transport(transport)

    print 'successful connection'

    return sftp, transport

def ftpDisconnect(ftp, trans):
    ftp.close()
    trans.close()

def ftpFile(subdir, file, ftp):
    print 'Ftp-ing: ' + str(file)

    path = '/external/avis/data/lockst/' 
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directy already exists...'
    path += subdir
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directory already exists...'
        
    ftp.put(os.getcwd() + '/' + file, path+'/'+file)


PASSWORD = ''

if __name__ == '__main__':
    
    global PASSWORD

    try:
        ASG = sys.argv[1]
        PASSWORD = str(sys.argv[2])
    except:
        print 'Usage: arg1 = path arg2 = ftp password'
        

    os.chdir(ASG)
    ls_output = Popen(['ls'], stdout=PIPE)

    #print ls_output.communicate()
    dirs = ls_output.communicate()[0].split('\n')

    ftp, trans = ftpConnect()

    for d in dirs:
        try:
            os.chdir(d)
        except:
            print 'Directory error for: ' + str(d)
            continue
        
        ls = Popen(['ls'], stdout=PIPE)
        vtp_list = ls.communicate()[0].split('\n')
        
        for vtp in vtp_list:
            try:
                 ftpFile(d,vtp, ftp)
            except:
                print 'sFTP error for: ' + str(vtp)
                continue
            
        os.chdir('..')

    ftpDisconnect(ftp,trans)

    print 'Done!'






