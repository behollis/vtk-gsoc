import memberReader
import numpy
from vtk.util import numpy_support as nsup
import vtk

MEMBERS = 200

def makeReadFunction(source, index):
    def readData():
        # Temporary hack. Need to make VTK hold reference to numpy
        global u, v, rho
        #print "executing with ", index
        u, v, rho = memberReader.readMember(index)
        npts = u.shape[0] * u.shape[1]

        vel = numpy.zeros((npts,3))
        vel[:,0] = u.reshape(npts,)
        vel[:,1] = v[0:127,0:127].reshape(npts,)

        velArray = nsup.numpy_to_vtk(vel)
        velArray.SetName("velocity")

        rhoArray = nsup.numpy_to_vtk(rho[0:127,0:127].reshape((npts,)))
        rhoArray.SetName("rho")

        image = source.GetCurrentOutputInformation().Get(vtk.vtkDataObject.DATA_OBJECT())

        image.SetDimensions(127, 127, 1)
        image.GetPointData().AddArray(velArray)
        image.GetPointData().AddArray(rhoArray)

    return readData


class ProgrammableFilterDataContainer(object):
    def __init__(self, f):
        self.Filter = f
        self.Count = 0
        self.MaxIterations = MEMBERS
        self.Streamlines = []

def mkRequestUpdateExtent(container):
    def requestUpdateExtent():
        self = container.Filter
        self.GetInputInformation(0, 0).Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), container.Count)
        #print 'here', container.Count
    return requestUpdateExtent

def mkRequestData(container):
    def requestData():
        self = container.Filter

        inpt = self.GetInputDataObject(0, 0)

        pt = vtk.vtkPointSource()
        pt.SetNumberOfPoints(1)
        pt.SetCenter(29.,29.,0.)
        pt.SetRadius(0.)

        st = vtk.vtkStreamTracer()
        st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, 'velocity')
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
            
            calcPCA(comp1_array, comp2_array)
            
            print 'writing'
            w = vtk.vtkPolyDataWriter()
            w.SetInputConnection(append.GetOutputPort())
            w.SetFileName("spaghetti.vtk")
            w.Write()

            return
        container.Count += 1
    return requestData

def calcPCA(xarray, yarray):
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

if __name__ == '__main__':
    
    r = vtk.vtkEnsembleSource()
    aColumn = vtk.vtkIntArray()
    aColumn.SetName("Ensemble Index")
    
    for mem in range(1,MEMBERS+1):
        ps = vtk.vtkProgrammableSource()
        ps.GetStructuredPointsOutput()
        ps.GetExecutive().SetOutputData(0, vtk.vtkImageData())
        ps.SetExecuteMethod(makeReadFunction(ps, mem))
        r.AddMember(ps)
        aColumn.InsertNextValue(mem)
    
    table = vtk.vtkTable()
    table.SetNumberOfRows(MEMBERS)
    table.GetRowData().AddArray(aColumn)
    r.SetMetaData(table)
    
    pf = vtk.vtkProgrammableFilter()
    pf.SetInputConnection(r.GetOutputPort())
    
    container = ProgrammableFilterDataContainer(pf)
    pf.SetRequestUpdateExtentMethod(mkRequestUpdateExtent(container))
    pf.SetExecuteMethod(mkRequestData(container))
    pf.Update()

