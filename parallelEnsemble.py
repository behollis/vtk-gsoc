import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys
import os
#from mpi4py import MPI
import h5py


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

NUM_CORES_PER_NODE = 5
NUM_NODES = 3
MEMBERS = 1000
SLICE = MEMBERS / NUM_CORES_PER_NODE #division of members per thread
EXT_X = 125
EXT_Y = 125
DATA_ROOT = '/home/behollis/DATA/out/'
            
if __name__ == '__main__':
    
    wc = vtk.vtkMPIController()
    
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
        
    if grank == 1:
        print 'terminating process #1'
        exit()
    
    grank = 0
    gsize = 1
    localSize = gsize / NUM_CORES_PER_NODE
    localGroup = 0#grank / localSize
    SLICE = 2; MEMBERS = 20
    
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
    
    #st.SetController(vtk.vtkDummyController())
    st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
    st.SetInputData(r.GetOutputDataObject(0))
   
    st.SetMaximumPropagation(500)
    st.SetMaximumNumberOfSteps(200) 
    st.SetInitialIntegrationStep(0.2)
    st.SetMinimumIntegrationStep(0.1)
    st.SetMaximumIntegrationStep(0.2)
    st.SetTerminalSpeed(0.001)
    st.SetIntegratorTypeToRungeKutta45()
    st.SetIntegrationDirection( vtk.vtkStreamTracer.BACKWARD )
    st.SetComputeVorticity(False)
    
    dir = DATA_ROOT+'lockExSt/ts00050/'
    if not os.path.exists(dir):
        os.makedirs(dir)
        
    #determine node
    node = grank / NUM_CORES_PER_NODE
    print 'node: ' + str(node)
    
    fname = str(node) + '.hdf5'
    f = h5py.File( dir+fname, 'w' )
        
    lrank = grank % NUM_CORES_PER_NODE
    
    #fill out grid of seeds for streamtracer
    colcells = EXT_Y 
    
    grid = vtk.vtkUnstructuredGrid()        
    pts = vtk.vtkPoints()
    
    #tot_pts = 0
    for x in range(node,int((node+1)*EXT_X/3.) + 1): 
        for y in range(node,(node+1)*EXT_Y + 1):
            print x
            print y
            #id = y + (y * x)
            
            #id = x + (y * rowcells)
            #cell_id = y + (x * colcells) 
            #print 'cell_id for pts ids: ' + str(cell_id)
            #pts.InsertPoint(cell_id, float(x),float(y),0.)
            #print 'pt id: ' + str(tot_pts)
            #tot_pts += 1
            pts.InsertNextPoint(float(x),float(y),0.)  
            #id += 1 
           
            
    grid.SetPoints(pts)
    st.SetSourceData(grid)  
    
    w = vtk.vtkPolyDataWriter()
    w.SetInputConnection( st.GetOutputPort() )
    
    #generate streamlines for each slice of the corresponding member vector field
    for mem in range( lrank*SLICE + 1, (lrank+1)*SLICE + 1 ):
        # set current member
        outInfo = r.GetOutputInformation(0)
        outInfo.Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), mem)
        r.Update()
        
        st.Update()
        
        w.Update()
        vtkPolyData_vtp = w.GetInput()
        
        #print 'polydata num cells: ' + str(vtkPolyData_vtp.GetNumberCells())
        
        print 'node %d' % node
        print 'EXT_X ' + str(int(EXT_X/3.))
        
        #cnt = vtkPolyData_vtp.GetNumberOfCells()
        #print 'number of cells: ' + str(cnt)
        
        for x in range(node,int((node+1)*EXT_X/3.)): 
            for y in range(node,(node+1)*EXT_Y):
        
                cell_id = y + (x * colcells) 
                
                print 'cell_id: ' + str(cell_id)
                
                try:
                    vtkPolyLine_sl = vtkPolyData_vtp.GetCell( cell_id )
                    
                    #if id > 40:
                    #    exit()
                   
                    #vtkPolyLine_sl = vtkPolyData_vtp.GetNextCell()
                
                
                    vtkPoints_pts = vtkPolyLine_sl.GetPoints()
                    
                    
                    
                    #copy streamline points to numpy array
                    num_pts = vtkPoints_pts.GetNumberOfPoints()
                    slnp = numpy.ndarray(shape=(2,num_pts))
                    for pt_id in range(0,num_pts):
                        
                        pt_tuple = vtkPoints_pts.GetPoint(pt_id)
                        slnp[0][pt_id] = pt_tuple[0]
                        slnp[1][pt_id] = pt_tuple[1]
                        
                        if pt_id == 0:
                            #print 'id: ' + str(id)
                            x = int(pt_tuple[0])
                            y = int(pt_tuple[1])
                            dir = '/mem' + str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3) 
                        
                        #print pt_tuple
                        
                    #id += 1
                        
                    #print dir
                    
                    f[dir] = slnp
                except:
                    #exit()
                    print 'missing cell id: ' + str(cell_id)
                    continue
                
                
                
