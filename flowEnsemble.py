import memberReader
import numpy
from vtk.util import numpy_support as nsup
import vtk


def makeReadFunction(source, index):
    def readData():
        # Temporary hack. Need to make VTK hold reference to numpy
        global u, v, rho
        print "executing with ", index
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
        self.MaxIterations = 2
        self.Streamlines = []

def mkRequestUpdateExtent(container):
    def requestUpdateExtent():
        self = container.Filter
        info = self.GetInputInformation(0, 0).Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), container.Count)
        print 'here', container.Count
    return requestUpdateExtent

def mkRequestData(container):
    def requestData():
        self = container.Filter

        input = self.GetInputDataObject(0, 0)

        line = vtk.vtkLineSource()
        line.SetPoint1(0, 0, 0)
        line.SetPoint2(126, 126, 0)
        line.SetResolution(10)

        st = vtk.vtkStreamTracer()
        st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
        st.SetSourceConnection(line.GetOutputPort())
        st.SetInputData(input)
        st.SetMaximumPropagation(100)
        st.SetInitialIntegrationStep(0.2)
        st.SetMinimumIntegrationStep(0.01)
        st.SetMaximumIntegrationStep(0.5)
        st.SetIntegratorTypeToRungeKutta45()
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
            for line in container.Streamlines:
                append.AddInputData(line)
            append.Update()

            print 'writing'
            w = vtk.vtkPolyDataWriter()
            w.SetInputConnection(append.GetOutputPort())
            w.SetFileName("spaghetti.vtk")
            w.Write()

            return
        container.Count += 1
    return requestData

if __name__ == '__main__':
    ps1 = vtk.vtkProgrammableSource()
    ps1.GetStructuredPointsOutput()
    ps1.GetExecutive().SetOutputData(0, vtk.vtkImageData())
    ps1.SetExecuteMethod(makeReadFunction(ps1, 1))
    
    ps2 = vtk.vtkProgrammableSource()
    ps2.GetStructuredPointsOutput()
    ps2.GetExecutive().SetOutputData(0, vtk.vtkImageData())
    ps2.SetExecuteMethod(makeReadFunction(ps2, 2))
    
    r = vtk.vtkEnsembleSource()
    
    aColumn = vtk.vtkIntArray()
    aColumn.SetName("Ensemble Index")
    for res in [1,2]:
        aColumn.InsertNextValue(res)
    table = vtk.vtkTable()
    table.SetNumberOfRows(2)
    table.GetRowData().AddArray(aColumn)
    r.SetMetaData(table)
    
    r.AddMember(ps1)
    r.AddMember(ps2)
    
    pf = vtk.vtkProgrammableFilter()
    pf.SetInputConnection(r.GetOutputPort())
    
    container = ProgrammableFilterDataContainer(pf)
    pf.SetRequestUpdateExtentMethod(mkRequestUpdateExtent(container))
    pf.SetExecuteMethod(mkRequestData(container))
    pf.Update()