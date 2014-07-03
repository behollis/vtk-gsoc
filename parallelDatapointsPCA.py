import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys

EXT_X = 40
EXT_Y = 40
MEMBERS = 20
DATAPOINTS = EXT_X * EXT_Y

# mpiexec -n <NUM_PROCESSES> pvtkpython parallelDatapointsPCA.py
# <NUM_PROCESSES> should be set to DATAPOINTS

def calcPCA(xarray, yarray):
    ''' Returns the eigenvalue the covariance 
        matrix for the terminal particle positions. 
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

def numpy2vtk(nparr, name):
    ''' Converts numpy array into a vtk double array '''
    out = vtk.vtkDoubleArray()
    out.SetNumberOfComponents(1)
    
    for ent in range(0,nparr.shape[0]):
        out.SetName( name )
        out.InsertNextValue(nparr[ent])
    
    return out

if __name__ == '__main__':
    
    wc = vtk.vtkMPIController()
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
    
    # Split the process space to 2 (because there are 2
    # ensemble members. Can be generalized to n members)
    # Note that if localSize > 1, there will be redundant
    # IO because of parallelization over seeds.
    localSize = gsize / DATAPOINTS
    localGroup = grank / localSize
    contr = None
    
    for i in range(DATAPOINTS):
        group = vtk.vtkProcessGroup()
        group.SetCommunicator(wc.GetCommunicator())
        for j in range(localSize):
            group.AddProcessId(i*localSize + j)
        c = wc.CreateSubController(group)
        if c != None:
            contr = c
    
    rank = contr.GetLocalProcessId()
    size = contr.GetNumberOfProcesses()
    
    x = grank % EXT_X
    y = grank / EXT_X
   
    fpts = numpy.ndarray(shape=(2,MEMBERS))
   
    for mem in range(0,MEMBERS):
        curr_x_mem = numpy.loadtxt('./out/x_member_%d.txt' % mem)
        curr_y_mem = numpy.loadtxt('./out/y_member_%d.txt' % mem)
        
        fpts[0][mem] = curr_x_mem[x][y]
        fpts[1][mem] = curr_y_mem[x][y]
    
    eigenvalue = calcPCA( numpy2vtk(fpts[0][:], 'x'), numpy2vtk(fpts[1][:], 'y') ) 
    
    f = open('./out/'+str(x)+'_'+str(y)+'_eigenvalue.txt', 'w+')
    f.write(str(eigenvalue))
    f.close() 
    
    
    
    
    
    