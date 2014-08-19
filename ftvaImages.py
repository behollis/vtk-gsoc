'''
Reads FTVA data into a vtkImageData object for viewing.
'''

import numpy
import vtk
import os
import sys
from subprocess import *

NUM_CORES = 1
INTERVAL = 20#intervals of steps where ftva is taken
STEPS = 500
SLICE = (STEPS / INTERVAL) / NUM_CORES #division of members per thread
X_EXT = 125
Y_EXT = 125
NULL_FTVA = -10
DATA_ROOT = '/home/nfs2/'

'''
wc = vtk.vtkMPIController()
gsize = wc.GetNumberOfProcesses()
grank = wc.GetLocalProcessId()
if gsize % 2 != 0:
    if grank == 0:
        print 'Only an even number of ranks is support'
    sys.exit(0)

localSize = gsize / NUM_CORES
localGroup = grank / localSize
'''

#print 'SLICE: ' + str(SLICE) + ' for grank: ' + str(grank)

grank = 0

dir = DATA_ROOT+'lockExSt/ts00050/ftva/'

vtkImage = vtk.vtkImageData()
vtkImage.SetExtent( 0, X_EXT, 0, Y_EXT, 0, SLICE)
#vtkImage.SetDimensions(X_EXT, Y_EXT, SLICE) 
vtkImage.AllocateScalars(vtk.VTK_DOUBLE, 1)

'''
vtkImageCells = vtk.vtkImageData()
vtkImageCells.SetExtent( 1, X_EXT, 1, Y_EXT, 1, SLICE )
#vtkImage.SetDimensions(X_EXT, Y_EXT, SLICE) 
vtkImageCells.AllocateScalars(vtk.VTK_DOUBLE, 1)
'''

vtkDArrayPts = vtk.vtkDoubleArray()
vtkDArrayPts.SetNumberOfComponents(1)
vtkDArrayPts.SetNumberOfTuples(vtkImage.GetNumberOfPoints())
vtkDArrayPts.SetName('ftvaPts')

vtkDArrayCells = vtk.vtkDoubleArray()
vtkDArrayCells.SetNumberOfComponents(1)
vtkDArrayCells.SetNumberOfTuples(vtkImage.GetNumberOfCells())
vtkDArrayCells.SetName('ftvaCells')
    
# fill block of ftva data for a range of integration steps
for istep in range( grank*SLICE, (grank+1)*SLICE + 1):
    for y in range(0,Y_EXT):
        for x in range(0,X_EXT):
            
            #print x
            #print y
            #print istep
            
            sdir = dir + 'x' + str(x).zfill(3) + '/'  
            ssdir = sdir + 'y' + str(y).zfill(3) + '/' 
           
            f = open(ssdir+'ftva.txt')
            lines = f.readlines()
            
            ftva = float(NULL_FTVA)
            
            if len(lines) >= istep+1:
                if ftva is not 'nan':
                    ftva = float(lines[istep])
                
            #print lines
            #print 'ftva: ' + str(ftva)
            #print x
            #print y
            #print istep
            '''
            if float(ftva) > 2.0:
                print istep
                print x
                print y
                exit()
            '''
            
            vtkDArrayPts.SetValue((istep+1)*(x+1)*(y+1), float(ftva))
            
            if x < (X_EXT - 1) and y < (Y_EXT - 1) and istep < ((grank+1)*SLICE - 1):
                vtkDArrayCells.SetValue((istep+1)*(x+1)*(y+1), float(ftva))
            
    print 'completed step: ' + str(istep)

vtkImage.GetPointData().AddArray(vtkDArrayPts)
vtkImage.GetCellData().AddArray(vtkDArrayCells) # for volumetrics

block_filename = DATA_ROOT + str(grank) + '_ftva.vti'

# Write block of image data to xml vtk format
writer = vtk.vtkXMLImageDataWriter()
writer.SetFileName(block_filename)
writer.SetInputData(vtkImage)
writer.Write()

print 'finished!'
    