import memberReader
import numpy
from vtk.util import numpy_support as nsup
from vtk.util import vtkAlgorithm as vta
import vtk
import sys

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
localSize = gsize / 2
localGroup = grank / localSize
contr = None

for i in range(2):
    group = vtk.vtkProcessGroup()
    group.SetCommunicator(wc.GetCommunicator())
    for j in range(localSize):
        group.AddProcessId(i*localSize + j)
    c = wc.CreateSubController(group)
    if c != None:
        contr = c

rank = contr.GetLocalProcessId()
size = contr.GetNumberOfProcesses()

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


ps1 = vtk.vtkPythonAlgorithm()
ps1.SetPythonObject(EnsembleReader(1))

ps2 = vtk.vtkPythonAlgorithm()
ps2.SetPythonObject(EnsembleReader(2))

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

r.UpdateInformation()
outInfo = r.GetOutputInformation(0)
outInfo.Set(vtk.vtkEnsembleSource.UPDATE_MEMBER(), localGroup)
r.Update()

nSeeds = 100
nLocalSeeds = nSeeds / size
start = rank*126.0/size
end = start + 126.0/size
line = vtk.vtkLineSource()
line.SetPoint1(start, start, 0)
line.SetPoint2(end, end, 0)
line.SetResolution(nLocalSeeds)

st = vtk.vtkStreamTracer()
# To disable the inherent data parallelism in the filter
st.SetController(vtk.vtkDummyController())
st.SetInputArrayToProcess(0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, "velocity")
st.SetSourceConnection(line.GetOutputPort())
st.SetInputData(r.GetOutputDataObject(0))
st.SetMaximumPropagation(100)
st.SetInitialIntegrationStep(0.2)
st.SetMinimumIntegrationStep(0.01)
st.SetMaximumIntegrationStep(0.5)
st.SetIntegratorTypeToRungeKutta45()

w = vtk.vtkXMLPolyDataWriter()
w.SetInputConnection(st.GetOutputPort())
w.SetFileName("slines%d.vtp" % grank)
w.Write()