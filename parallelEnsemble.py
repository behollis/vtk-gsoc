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

# mpiexec -n <NUM_CORES> pvtkpython parallelEnsemble.py

NUM_CORES = 12
MEMBERS = NUM_CORES * 83
SLICE = MEMBERS / NUM_CORES #division of members per thread
EXT_X = 127
EXT_Y = 127
ROOT = '/home/behollis/'

if __name__ == '__main__':
    
    wc = vtk.vtkMPIController()
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
   
    localSize = gsize / NUM_CORES
    localGroup = grank / localSize
    
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
    
    '''
    grid = vtk.vtkUnstructuredGrid()
    pts = vtk.vtkPoints()
    pts.InsertPoint(float(x),float(y),0.)   
    grid.SetPoints(pts)
    '''
   
    #st.SetSourceConnection( pt.GetOutputPort() )
    st.SetMaximumPropagation(5000)
    st.SetMaximumNumberOfSteps(1000) 
    st.SetInitialIntegrationStep(0.2)
    st.SetMinimumIntegrationStep(0.1)
    st.SetMaximumIntegrationStep(0.2)
    st.SetTerminalSpeed(0.000001)
    st.SetIntegratorTypeToRungeKutta45()
    st.SetIntegrationDirectionToBoth()
    
    dir = ROOT+'lockExSt/ts00050/'
    if not os.path.exists(dir):
        os.makedirs(dir)
        
    for mem in range( grank*SLICE, (grank+1)*SLICE ):
        print 'mem: ' + str(mem).zfill(4)
        for x in range(0,EXT_X):
            
            sdir = dir + str(mem).zfill(4) + '/x' + str(x).zfill(3) 
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
