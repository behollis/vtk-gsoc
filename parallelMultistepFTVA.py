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
    
    '''
    #forward sl
    vtkPolyLine_sl = vtkPolyData_vtp.GetCell(0)
    vtkPoints_pts = vtkPolyLine_sl.GetPoints()
    
    stlst_f.append(vtkPoints_pts)
    '''
    
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
        
        
        try:
            ls = Popen(['ls'], stdout=PIPE)
            vtps = list()
            for f in ls.communicate()[0].split('\n'):
                if '.vtp' in f:
                    vtps.append(f)
    
            #print 'reading streamline: ' + str(mem_dir)
            if len(vtps) is not 0:
                readStreamline(vtps[0], stlst_f, stlst_b)
        except:
            print 'Error reading streamline: ' + str(mem_dir)
                 
        os.chdir('../../..') # go to next member
    
def performStats(x,y,st_f, st_b):
    #print 'performing stats on: ' + str(x) + ', ' + str(y)
    
    dir = OUT_PATH + 'ftva/x' + str(x).zfill(3) + '/y' + str(y).zfill(3) 
    if not os.path.exists(dir):
        os.makedirs(dir)
    
    os.chdir(dir)
    
    ls_output = Popen(['ls'], stdout=PIPE)
    files = ls_output.communicate()[0].split('\n')
    
    out = None
    #if 'ftva' in files:
    #    out = open('ftva', 'a')
    #else:
    out = open('ftva.txt', 'w+')
    
    for pt_id in range(0,1000, INTERVAL):
        comp1_array = vtk.vtkDoubleArray()
        comp1_array.SetNumberOfComponents(1)
        comp2_array = vtk.vtkDoubleArray()
        comp2_array.SetNumberOfComponents(1)
        
        for line in st_b:
            pt_tuple = line.GetPoint(pt_id)
            
            comp1_array.SetName( 'x' )
            comp1_array.InsertNextValue(pt_tuple[0])
         
            comp2_array.SetName( 'y' )
            comp2_array.InsertNextValue(pt_tuple[1])
                
        var = calcPCA(comp1_array, comp2_array)
        
        #print 'writing ftva @ step: ' + str(pt_id) + ' for: ' + str(x) + ' ' + str(y)
        out.write(str(var)+'\n')
        
    out.close()
    
    
def main():
    for x in range(X_STR, X_END, SEED_RES):
        for y in range(Y_STR, Y_END, SEED_RES):
            print 'collecting streamlines for: ' + str(x) + ' ' + str(y)
            collectEnsembleStreamlines(x, y, SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            try:
                performStats(x,y,SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            except:
                print 'Error computing PCA for: ' + str(x) + ' ' + str(y)
            del SEED_STREAMLINES_F[:]; del SEED_STREAMLINES_B[:]
            
            
if __name__ == '__main__':
    PATH = '/home/behollis/lockExSt/ts00050/'
    OUT_PATH = '/media/behollis/TOSHIBA EXT/'
    SEED_STREAMLINES_F = list()
    SEED_STREAMLINES_B = list()
    X_EXT = 125
    Y_EXT = 125
    NUM_CORES = 4
    SLICE = X_EXT / 2 #X_EXT / NUM_CORES #division of members per thread
    
    SEED_RES = 1 # regularity of seed sampling for FTVA
    INTERVAL = 1 # interval of integration steps for multi-step FTVA
    
    wc = vtk.vtkMPIController()
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
    
    localSize = gsize / NUM_CORES
    localGroup = grank / localSize
    
    # using four cores
    B_OFFSET = 10
    E_OFFSET = 10
    X_STR = 10; Y_STR = 10
    X_END = X_EXT - 10; Y_END = Y_EXT - 10
    
    if grank == 0:
        X_STR = B_OFFSET
        X_END = X_EXT / 2
        Y_STR = B_OFFSET
        Y_END = Y_EXT / 2
    elif grank == 1:
        X_STR = X_EXT / 2
        X_END = X_EXT - E_OFFSET
        Y_STR = B_OFFSET
        Y_END = Y_EXT / 2
    elif grank == 2:
        X_STR = B_OFFSET
        X_END = X_EXT / 2
        Y_STR = Y_EXT / 2
        Y_END = Y_EXT - E_OFFSET
    else: #grank == 3:
        X_STR = X_EXT / 2
        X_END = X_EXT - E_OFFSET
        Y_STR = Y_EXT / 2
        Y_END = Y_EXT - E_OFFSET
        
    main()
    
    
    
    
    
    

    
    
    
       
