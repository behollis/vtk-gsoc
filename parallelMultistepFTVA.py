import numpy
import vtk
import os
import sys
from subprocess import *

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
    
    return eigenvalues.GetValue(0)


def readStreamline(stfile,stlst_f, stlst_b):
    
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(stfile)
    reader.Update()

    #one streamline per vtp
    vtkPolyData_vtp = reader.GetOutput()
    
    #forward sl
    vtkPolyLine_sl = vtkPolyData_vtp.GetCell(0)
    vtkPoints_pts = vtkPolyLine_sl.GetPoints()
    
    stlst_f.append(vtkPoints_pts)
    
    #backward sl
    vtkPolyLine_sl = vtkPolyData_vtp.GetCell(1)
    vtkPoints_pts = vtkPolyLine_sl.GetPoints()
    
    stlst_b.append(vtkPoints_pts)
   

def collectEnsembleStreamlines(x, y, stlst_f, stlst_b):
    ''' Read all points for each streamline thru a seed for the ensemble. '''
    
    os.chdir(PATH)
    ls_output = Popen(['ls'], stdout=PIPE)

    #print ls_output.communicate()
    dirs = ls_output.communicate()[0].split('\n')
    
    for mem_dir in dirs:
        try:
            os.chdir( mem_dir + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3) )
        except:
            print 'Directory error for: ' + str(mem_dir)
            continue
        
        
        ls = Popen(['ls'], stdout=PIPE)
        vtps = list()
        for f in ls.communicate()[0].split('\n'):
            if '.vtp' in f:
                vtps.append(f)
        
        readStreamline(vtps[0], stlst_f, stlst_b)
                 
        os.chdir('../../..') # go to next member
    
def performStats(x,y):
    print 'performing stats on: ' + str(x) + ', ' + str(y)
    return
    
def main():
    for x in range(40,X_EXT):
        for y in range(40, Y_EXT):
            collectEnsembleStreamlines(x, y, SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            performStats(x,y,SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            del SEED_STREAMLINES_F[:]; del SEED_STREAMLINES_B[:]
            
if __name__ == '__main__':
    PATH = '/home/behollis/lockExSt/ts00050/'
    SEED_STREAMLINES_F = list()
    SEED_STREAMLINES_B = list()
    X_EXT = 127
    Y_EXT = 127
    
    main()
    
    
    
    
    
    

    
    
    
       
