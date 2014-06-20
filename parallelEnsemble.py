import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys


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
    
class PCAFilter(vta.VTKAlgorithm):
    def __init__(self):
        vta.VTKAlgorithm.__init__(self, nInputPorts=1, outputType='vtkImageData')
        self.Count = 0
        self.MaxIterations = MEMBERS
        self.Streamlines = []

    def RequestInformation(self, vtkself, request, inInfo, outInfo):
        vtkSDDP = vtk.vtkStreamingDemandDrivenPipeline
        outInfo.GetInformationObject(0).Set(vtkSDDP.WHOLE_EXTENT(), (0, 126, 0, 126, 0, 0), 6)
        return 1
    
    def CalcPCA(self, xarray, yarray):
        ''' Computes the eigenvalues/eigenvectors of the covariance matrix for the terminal particle positions. '''
        
        #See: http://www.vtk.org/Wiki/VTK/Examples/Cxx/Utilities/PCAStatistics
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
     
        # Eigenvalues 
        eigenvalues = vtk.vtkDoubleArray()
        pcaStatistics.GetEigenvalues(eigenvalues)
    
        for idx in range(0,eigenvalues.GetNumberOfTuples()):
        #for eigenvalue in eigenvalues:
            print 'Eigenvalue ' + str(idx) + ' = ' + str(eigenvalues.GetValue(idx))
     
        # Eigenvectors 
        eigenvectors = vtk.vtkDoubleArray()
        pcaStatistics.GetEigenvectors(eigenvectors)
        
        for idx in range(0,eigenvectors.GetNumberOfTuples()):
            print 'Eigenvector ' + str(idx) + ' : '
            evec = [0]*eigenvectors.GetNumberOfComponents()
            eigenvectors.GetTuple(idx, evec)
            print evec

    def RequestData(self, vtkself, request, inInfo, outInfo):
        contr = vtk.vtkMPIController()
        rank = contr.GetLocalProcessId()

        inpt = self.GetInputDataObject(0, 0)

        pt = vtk.vtkPointSource()
        pt.SetNumberOfPoints(1)
        pt.SetCenter(29.,29.,0.)
        pt.SetRadius(0.)

        st = vtk.vtkPStreamTracer()
        
        st.SetController(vtk.vtkDummyController())
        st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
        st.SetSourceConnection(pt.GetOutputPort())
        st.SetInputData(inpt)
        
        #Specify the maximum length of a streamline expressed in LENGTH_UNIT. 
        st.SetMaximumPropagation(5000)
        
        st.SetMaximumNumberOfSteps(100) #bifurcation @ 100 steps
        st.SetInitialIntegrationStep(0.5)
        st.SetMinimumIntegrationStep(0.5)
        st.SetMaximumIntegrationStep(1.0)
        st.SetIntegratorTypeToRungeKutta45()
        st.SetIntegrationDirectionToBackward()#Forward()
        st.Update()

        sline = st.GetOutput().NewInstance()
        sline.ShallowCopy(st.GetOutput())

        container.Streamlines.append(sline)

        req = self.GetCurrentRequest()
        if container.Count == 0:
            req.Set(vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING(), 1)
        elif container.Count == container.MaxIterations - 1:
            container.Count = 0
            req.Remove(vtk.vtkStreamingDemandDrivenPipeline.CONTINUE_EXECUTING())
            append = vtk.vtkAppendPolyData()
            
            comp1_array = vtk.vtkDoubleArray()
            comp1_array.SetNumberOfComponents(1)
            comp2_array = vtk.vtkDoubleArray()
            comp2_array.SetNumberOfComponents(1)
            
            for line in container.Streamlines:
                append.AddInputData(line)
                
                tpt = line.GetPoint(line.GetNumberOfPoints()-1)
                
                comp1_array.SetName( 'x' )
                comp1_array.InsertNextValue(tpt[0])
                
                comp2_array.SetName( 'y' )
                comp2_array.InsertNextValue(tpt[1])
               
            append.Update()
            
            #CalcPCA(comp1_array, comp2_array)
            
            w = vtk.vtkXMLPolyDataWriter()
    
            w.SetInputConnection(st.GetOutputPort())
            w.SetFileName("slines%d.vtp" % grank)
            w.Write()

            return
        
        container.Count += 1

# mpiexec -n <NUM_PROCESSES> pvtkpython parallelEnsemble.py
# this should be less than or equal to the number of processes
MEMBERS = 20

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
    localSize = gsize / MEMBERS
    localGroup = grank / localSize
    contr = None
    
    for i in range(MEMBERS):
        group = vtk.vtkProcessGroup()
        group.SetCommunicator(wc.GetCommunicator())
        for j in range(localSize):
            group.AddProcessId(i*localSize + j)
        c = wc.CreateSubController(group)
        if c != None:
            contr = c
    
    rank = contr.GetLocalProcessId()
    size = contr.GetNumberOfProcesses()
    
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
    
    # Add programmable filter as vtkPythonAlgorithm
    pf = vtk.vtkPythonAlgorithm()
    
    pf.SetPythonObject(PCAFilter())
    pf.SetInputConnection(r.GetOutputPort())
    pf.Update()
    
   