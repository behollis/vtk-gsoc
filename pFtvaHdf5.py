import vtk
import os
import sys
import time
import h5py
import numpy as np
from vtk.util import numpy_support as nsup

from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank

DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'

'''
Parallel FTVA using MPI and parallel hdf5.
Started: October 13th, 2014 (Monday).
Author: Brad Hollister 
'''

def calcPCA(xarray, yarray):
    ''' Returns the eigenvalue of the covariance matrix for the terminal particle 
        positions. 
    '''
    datasetTable = vtk.vtkTable()
    datasetTable.AddColumn(xarray)
    datasetTable.AddColumn(yarray)
    
    pcaStatistics = vtk.vtkPCAStatistics()
    
    pcaStatistics.SetInputData( vtk.vtkStatisticsAlgorithm.INPUT_DATA, datasetTable )
    
    pcaStatistics.SetColumnStatus('x', 1 )
    pcaStatistics.SetColumnStatus('y', 1 )
    pcaStatistics.RequestSelectedColumns()
    pcaStatistics.SetDeriveOption(True)
    pcaStatistics.Update()
     
    eigenvalues = vtk.vtkDoubleArray()
    pcaStatistics.GetEigenvalues(eigenvalues)
    
    return eigenvalues
  
def main(dset):
    for x in range( node*(X_EXT/3), (node+1)*(X_EXT/3)+1 ): #range( 0, X_EXT ):
        if rank == 1:
            print 'x: ' + str(x)
        for y in range( rank*SLICE_Y, (rank+1)*SLICE_Y ):
                        
            xarray = vtk.vtkFloatArray()
            xarray.SetName( 'x' )
            yarray = vtk.vtkFloatArray()
            yarray.SetName( 'y' )
            
            for step in range(0,STEPS,INTERVAL):
                for mem in range( 1, 150):#MEMBERS, MINTERVAL):
                    dir = '/mem' + str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3)
                  
                    #print 'step: ' + str(step)
                    #print 'mem: ' + str(mem)
                    
                    try:
                        pt0_0 = f0[dir][0][step]
                        pt1_0 = f0[dir][1][step]
                        
                        #print pt0_0
                        #print pt1_0
                    except:
                        #print 'could not read values from hdf5 file'
                        #print '\t for x={0} y={1} mem={2}'.format( x, y, mem )
                        continue
                    
                    xarray.InsertNextValue(pt0_0)
                    yarray.InsertNextValue(pt1_0)
                
                    #print 'r={0} x={1} y={2} mem={3} step={4}'.format(rank, x, y, mem, step)
                    #print 'first  streamline point: ' + str(pt0_0)
                    #print 'second streamline point: ' + str(pt1_0)
                    
                #get PCA for this step...
                evals = calcPCA(xarray, yarray)
                majeval = 0.0
                try:
                    majeval = evals.GetValue(0)
                except:
                    print 'missing majeval'
                    
                if step % INTERVAL == 0:
                    dset[x,y,step/INTERVAL] = majeval
                
       
NUM_PROC = 6 #number of processes / cores per node
NUM_NODES = 1
X_EXT = 125
Y_EXT = 125 
SLICE_Y = int( Y_EXT / float( NUM_PROC ) )
SLICE_X = int ( X_EXT / float( NUM_NODES ) )
STEPS = 300
INTERVAL = 10
MEMBERS = 1000
MINTERVAL = 100
        
if __name__ == '__main__':
    
    comm = MPI.COMM_WORLD
    rank = comm.rank
    #rank = 0
    
    #determine node
    node = 0#rank / NUM_PROC
    print 'node: ' + str(node)
    
    #print 'reading hdf5 file...'
    f0 = h5py.File(DPATH+'0.0.hdf5', 'r', driver='mpio', comm=comm)
    #f1 = h5py.File(DPATH+'0.1.hdf5', driver='mpio', comm=comm)
    #f2 = h5py.File(DPATH+'0.2.hdf5', driver='mpio', comm=comm)
    #f3 = h5py.File(DPATH+'0.3.hdf5', driver='mpio', comm=comm)
    #f4 = h5py.File(DPATH+'0.4.hdf5', driver='mpio', comm=comm)
    #f5 = h5py.File(DPATH+'0.5.hdf5', driver='mpio', comm=comm)
    #print 'finished reading h5py file!'
    
    output = h5py.File( DPATH+'ftva.hdf5', 'w', driver='mpio', comm=comm)
    
    dset = output.create_dataset(name='majEigenVal', \
                                 shape=(X_EXT/3 + 1,Y_EXT,(STEPS/INTERVAL)+1), dtype='f')
    
                
    main(dset)
    f0.close()
    output.close()
    
    print 'rank: ' + str(rank) + ' finished!'
    
    
    
    
   
    