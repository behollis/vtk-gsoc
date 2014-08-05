'''
Reads FTVA data into a vtkImageData object for viewing.
'''

import numpy
import vtk
import os
import sys
from subprocess import *

NUM_CORES = 25
STEPS = 1000
SLICE = STEPS / NUM_CORES #division of members per thread
EXT_X = 125
EXT_Y = 125
DATA_ROOT = '/home/nfs2/'

wc = vtk.vtkMPIController()
gsize = wc.GetNumberOfProcesses()
grank = wc.GetLocalProcessId()
if gsize % 2 != 0:
    if grank == 0:
        print 'Only an even number of ranks is support'
    sys.exit(0)

localSize = gsize / NUM_CORES
localGroup = grank / localSize

dir = DATA_ROOT+'lockExSt/ts00050/'

vtkImage = vtk.vtkImageData()
vtkImage.SetDimensions(X_EXT, Y_EXT, SLICE)
    
# fill block of ftva data for a range of integration steps
for istep in range( grank*SLICE, (grank+1)*SLICE ):
    
    for x in range(0,EXT_X):
        
        sdir = dir + '/x' + str(x).zfill(3) 
        if not os.path.exists(sdir):
            os.makedirs(sdir)
        
        os.chdir(sdir)
        
        for y in range(0,EXT_Y):
            
            ssdir = sdir + '/y' + str(y).zfill(3) 
            
            os.chdir(ssdir)
            
            f = open('ftva.txt')
            lines = f.readlines()
            
            ftva = lines[istep]
            
            #Example usage: http://www.vtk.org/Wiki/VTK/Examples/Cxx/IO/WriteVTI
            pixel = vtkImage.GetScalarPointer(x,y,istep)
            pixel[0] = ftva
            
            os.chdir('..')

# Write block of image data to xml vtk format
writer = vtk.vtkXMLImageDataWriter()
block_filename = str(grank) + '_ftva.vti'
writer.SetFileName(dir+block_filename)
writer.SetInputData(vtkImage)
writer.Write()

print 'finished!'
    