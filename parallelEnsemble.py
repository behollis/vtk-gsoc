import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys
import os
from mpi4py import MPI


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
        print 'reading member: ' + str(self.Index)
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

NUM_CORES = 4#20
MEMBERS = NUM_CORES * 5
SLICE = MEMBERS / NUM_CORES #division of members per thread
EXT_X = 127
EXT_Y = 127
DATA_ROOT = '/home/data_local/'

#c = vtk.vtkMultiProcessController.GetGlobalController()
#comm = vtk.vtkMPI4PyCommunicator.ConvertToPython(c.GetCommunicator())



def initRemainingMems():
    ''' Init unprocessed members for seeds counter.
    '''
    seeds = dict()
    for m in range(0,MEMBERS):
        for x in range(0,2):#EXT_X):
            for y in range(0,2):#EXT_Y):
                #unprocessed members
                seeds[(x,y)] = MEMBERS
                
    return seeds
            
if __name__ == '__main__':
    '''
    wc = vtk.vtkMPIController()
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
    '''
    
    comm = MPI.COMM_WORLD
    grank = comm.Get_rank()
   
    #localSize = gsize / NUM_CORES
    #localGroup = grank / localSize
    
    #http://mpi4py.scipy.org/docs/usrman/tutorial.html#collective-communication
    
    if grank == 0:
        #wc.Initialize()
        mem_seeds = initRemainingMems()
        mem_seeds = comm.bcast(mem_seeds, root=0)
    else:
        mem_seeds = None
    
    print 'grank: ' + str(grank)
    print mem_seeds
    
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
   
    st = vtk.vtkStreamTracer()
    st.SetController(vtk.vtkDummyController())
    st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
    st.SetInputData(r.GetOutputDataObject(0))
   
    st.SetMaximumPropagation(500)
    st.SetMaximumNumberOfSteps(200) 
    st.SetInitialIntegrationStep(0.2)
    st.SetMinimumIntegrationStep(0.1)
    st.SetMaximumIntegrationStep(0.2)
    st.SetTerminalSpeed(0.001)
    st.SetIntegratorTypeToRungeKutta45()
    st.SetIntegrationDirection(BACKWARD)
    st.SetComputeVorticity(false)
    
    dir = DATA_ROOT+'lockExSt/ts00050/'
    if not os.path.exists(dir):
        os.makedirs(dir)
        
    for mem in range( grank*SLICE + 1, (grank+1)*SLICE + 1 ):
        
        # set current member
        outInfo = r.GetOutputInformation(0)
        outInfo.Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), mem)
        r.Update()
        
        print 'mem: ' + str(mem).zfill(4)
        for x in range(0,EXT_X):
            
            sdir = dir + 'mem' + str(mem).zfill(4) + '/x' + str(x).zfill(3) 
            if not os.path.exists(sdir):
                os.makedirs(sdir)
            
            os.chdir(sdir)
            
            for y in range(0,EXT_Y):
                
                ssdir = sdir + '/y' + str(y).zfill(3) 
                if not os.path.exists(ssdir):
                    os.makedirs(ssdir)
            
                os.chdir(ssdir)
    
                w = vtk.vtkXMLPolyDataWriter()
                grid = vtk.vtkUnstructuredGrid()
                pts = vtk.vtkPoints()
            
                pts.InsertNextPoint(float(x),float(y),0.)   
                grid.SetPoints(pts)
                
                w.SetFileName('sline_M' + str(mem).zfill(4) + '_X' + str(x).zfill(3) + '_Y' + str(y).zfill(3) +'.vtp')
                st.SetSourceData(grid)  
                st.Update()
    
                w.SetInputConnection(st.GetOutputPort())
                w.Update()
                w.Write()
                
                os.chdir('..')
                
    '''
    #wc.Finalize()
