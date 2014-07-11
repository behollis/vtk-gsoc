import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys
import os


class EnsembleReader(vta.VTKAlgorithm):
    def __init__(self, index):
        vta.VTKAlgorithm.__init__(self, nInputPorts=0, outputType='vtkImageData')
        self.Index = index

    def RequestInformation(self, vtkself, request, inInfo, outInfo):
        vtkSDDP = vtk.vtkStreamingDemandDrivenPipeline
        outInfo.GetInformationObject(0).Set(vtkSDDP.WHOLE_EXTENT(), (0, 126, 0, 126, 0, 0), 6)
        return 1

    def RequestData(self, vtkself, request, inInfo, outInfo):
        contr = vtk.vtkMPIController()
        rank = contr.GetLocalProcessId()
        u, v, self.Rho = memberReader.readMember(self.Index)
        npts = u.shape[0] * u.shape[1]

        self.Vel = numpy.zeros((npts,3))
        self.Vel[:,0] = u.reshape(npts,)
        self.Vel[:,1] = v[0:127,0:127].reshape(npts,)

        velArray = nsup.numpy_to_vtk(self.Vel)
        velArray.SetName("velocity")

        rhoArray = nsup.numpy_to_vtk(self.Rho[0:127,0:127].reshape((npts,)))
        rhoArray.SetName("rho")

        image = outInfo.GetInformationObject(0).Get(vtk.vtkDataObject.DATA_OBJECT())

        image.SetDimensions(127, 127, 1)
        image.GetPointData().AddArray(velArray)
        image.GetPointData().AddArray(rhoArray)
        
        return 1

# mpiexec -n <NUM_PROCESSES> pvtkpython parallelEnsemble.py

NUM_CORES = 12#16#8
NUM_PROCS = 1
MEMBERS = NUM_CORES * 83#62#125
SLICE = MEMBERS / (NUM_CORES * NUM_PROCS) #division of members per thread
EXT_X = 127
EXT_Y = 127
ROOT = '../../../'

def calcPCA(xarray, yarray):
    ''' Returns the eigenvalue the covariance matrix for the terminal particle 
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

    #print 'gsize ' + str(gsize)
    #print 'grank ' + str(grank)
    localSize = gsize / NUM_CORES
    
    #print 'localsize ' + str(localSize)
    localGroup = grank / localSize
    
    contr = None
    
    '''
    for i in range(MEMBERS):
        group = vtk.vtkProcessGroup()
        group.SetCommunicator(wc.GetCommunicator())
        for j in range(localSize):
            group.AddProcessId(i*localSize + j)
        c = wc.CreateSubController(group)
        if c != None:
            contr = c
    
    rank = contr.GetLocalProcessId()
    print 'rank: ' + str(rank)
    size = contr.GetNumberOfProcesses()
    '''
    
    r = vtk.vtkEnsembleSource()
    
    table = vtk.vtkTable()
    table.SetNumberOfRows(MEMBERS)
    
    aColumn = vtk.vtkIntArray()
    aColumn.SetName("Ensemble Index")
   
    
    for mem in range(1,MEMBERS+1):
        # Add source / reader
        ps = vtk.vtkPythonAlgorithm()
        ps.SetPythonObject(EnsembleReader(mem))
        r.AddMember(ps)
         
        aColumn.InsertNextValue(mem)
        
    table.GetRowData().AddArray(aColumn)
    r.SetMetaData(table)
    r.UpdateInformation()
    outInfo = r.GetOutputInformation(0)
    outInfo.Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), localGroup)
    r.Update()
    
    #terninal points over field in a member
    #fpts_x = numpy.zeros(shape=(EXT_X, EXT_Y))
    #fpts_y = numpy.zeros(shape=(EXT_X, EXT_Y))
    
    st = vtk.vtkStreamTracer()
    
    st.SetController(vtk.vtkDummyController())
    st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
    st.SetInputData(r.GetOutputDataObject(0))
    
    # seed points for half of grid
    grid1 = vtk.vtkUnstructuredGrid()
    pts1 = vtk.vtkPoints()
    
    for x in range( 0, EXT_X-8, 4 ):
        for y in range( 0, EXT_Y-8, 4 ):
            pts1.InsertNextPoint(float(x),float(y),0.)
    
    # seed points for other half of grid
    grid2 = vtk.vtkUnstructuredGrid()
    pts2 = vtk.vtkPoints()
    
    for x in range( 3, EXT_X-8, 4 ):
        for y in range( 1, EXT_Y-8, 4 ):
            pts2.InsertNextPoint(float(x),float(y),0.)
            
   
    grid1.SetPoints(pts1)
    grid2.SetPoints(pts2)
   
    #st.SetSourceConnection( pt.GetOutputPort() )
    st.SetMaximumPropagation(5000)
    st.SetMaximumNumberOfSteps(500) 
    st.SetInitialIntegrationStep(0.5)
    st.SetMinimumIntegrationStep(0.1)
    st.SetMaximumIntegrationStep(0.5)
    st.SetTerminalSpeed(0.000001)
    st.SetIntegratorTypeToRungeKutta45()
    
    dir = ROOT+'lockExchangeStreamlinesTs0050/'+str(grank*SLICE)+'/'
    if not os.path.exists(dir):
        os.makedirs(dir)
   
    for mem in range( grank*SLICE, (grank+1)*SLICE ):
        # write out streamlines for this member(s)
        w = vtk.vtkXMLPolyDataWriter()
        
        mem_num = str(mem).zfill(4)
        print 'mem: ' + mem_num
        
        #1F
        w.SetFileName(dir+'1Fslines%s.vtp' % mem_num)
        st.SetSourceData(grid1)  
        st.SetIntegrationDirectionToForward()
        st.Update()
        
        w.SetInputConnection(st.GetOutputPort())
        w.Update()
        w.Write()
        
        #1B
        w.SetFileName(dir+'1Bslines%s.vtp' % mem_num)
        st.SetSourceData(grid1)  
        st.SetIntegrationDirectionToBackward()
        st.Update()
        
        w.SetInputConnection(st.GetOutputPort())
        w.Update()
        w.Write()
        
        #2F
        w.SetFileName(dir+'2Fslines%s.vtp' % mem_num)
        st.SetSourceData(grid2)
        st.SetIntegrationDirectionToForward()
        st.Update()
        
        w.SetInputConnection(st.GetOutputPort())
        w.Update()
        w.Write()
        
        #2B
        w.SetFileName(dir+'2Bslines%s.vtp' % mem_num)
        st.SetSourceData(grid2)
        st.SetIntegrationDirectionToBackward()
        st.Update()
        
        w.SetInputConnection(st.GetOutputPort())
        w.Update()
        w.Write()
        
        
       
        
        
        
        '''    
        line = st.GetOutput().NewInstance()
        line.ShallowCopy(st.GetOutput())
    
        tpt = line.GetPoint(line.GetNumberOfPoints()-1)
        
        fpts_x[x][y] = tpt[0]
        fpts_y[x][y] = tpt[1]
        '''
        
        #numpy.savetxt('./out/x_member_%d.txt' % grank, fpts_x)
        #numpy.savetxt('./out/y_member_%d.txt' % grank, fpts_y)
        
       
