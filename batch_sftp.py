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
    
    host ='moondance.soe.ucsc.edu'
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

def ftpFile(subdir,subsubxdir, subsubydir, file, ftp):
    print 'Ftp-ing: ' + str(file)

    #path = '/external/avis/data/lockst2/' 
    path = '.html/'
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directy already exists...'
        
    path += 'ts00050/' 
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directy already exists...'
        
    path += subdir
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directory already exists...'
        
    path += subsubxdir
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directory already exists...'
        
    path += subsubydir
    try:
        ftp.mkdir(path, mode=511)
    except:
        print 'directory already exists...'
        
    ftp.put(os.getcwd() + '/' + file, path+'/'+file)
    

PASSWORD = ''

if __name__ == '__main__':
    
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
    
    #dirs.reverse()

    for mem_dir in dirs:
        
        #print mem_dir
        
        #if '734' in mem_dir:
        #    print 'exiting on 734...'
        #    exit()
        
        try:
            os.chdir(mem_dir)
        except:
            print 'Directory error for: ' + str(mem_dir)
            continue
        
        #enter streamline dir
        ls0 = Popen(['ls'], stdout=PIPE)
        seed_list_x = ls0.communicate()[0].split('\n')
        
        for seed_dir_x in seed_list_x:
            
            if 'x027' not in seed_dir_x and 'x029' not in seed_dir_x:
                continue
            
            print 'got an x'
            
            try:
                os.chdir(seed_dir_x)
            except:
                print 'Directory error for: ' + str(seed_dir_x)
                continue
        
            ls1 = Popen(['ls'], stdout=PIPE)
            seed_list_y = ls1.communicate()[0].split('\n')
            
            for seed_dir_y in seed_list_y:
                
                if 'y028' not in seed_dir_y and 'y030' not in seed_dir_y:
                    continue
                
                print 'got a y'
            
                try:
                    os.chdir(seed_dir_y)
                except:
                    print 'Directory error for: ' + str(seed_dir_y)
                    continue
            
                ls2 = Popen(['ls'], stdout=PIPE)
                vtp_list = ls2.communicate()[0].split('\n')
                
                for vtp in vtp_list:
                    try:
                         ftpFile(mem_dir+'/', seed_dir_x+'/', seed_dir_y+'/', vtp, ftp)
                    except:
                        print 'sFTP error for: ' + str(vtp)
                        continue
                    
                os.chdir('..') #go to next y dir
                
            os.chdir('..') # go to next x dir
            
        os.chdir('..') # go to next member

    ftpDisconnect(ftp,trans)

    print 'Done!'






