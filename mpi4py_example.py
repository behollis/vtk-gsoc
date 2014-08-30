'''

import vtk

c = vtk.vtkMultiProcessController.GetGlobalController()

rank = c.GetLocalProcessId()

if rank == 0:
    ssource = vtk.vtkSphereSource()
    ssource.Update()
    c.Send(ssource.GetOutput(), 1, 1234)
else:
    sphere = vtk.vtkPolyData()
    c.Receive(sphere, 0, 1234)
    print sphere

'''

#http://www.kitware.com/blog/home/post/721
'''
import vtk

c = vtk.vtkMultiProcessController.GetGlobalController()

rank = c.GetLocalProcessId()

if rank == 0:
    ssource = vtk.vtkSphereSource()
    ssource.Update()
    sphere = ssource.GetOutput()
else:
    sphere = vtk.vtkPolyData()

c.Broadcast(sphere, 0)

if rank == 1:
    print sphere
    
    
from vtk.numpy_interface import dataset_adapter as dsa
from mpi4py import MPI
comm = vtk.vtkMPI4PyCommunicator.ConvertToPython(c.GetCommunicator())
import numpy

if rank == 0:
    ssource = vtk.vtkSphereSource()
    ssource.Update()
    ca = vtk.vtkCharArray()
    vtk.vtkCommunicator.MarshalDataObject(ssource.GetOutput(), ca)
    a = dsa.vtkDataArrayToVTKArray(ca)
    comm.send(a.shape[0], 1)
    comm.Send([a, MPI.CHAR], 1)
else:
    sz = comm.recv()
    buff = numpy.empty(sz, dtype=numpy.int8)
    comm.Recv([buff, MPI.CHAR])
    pd = vtk.vtkPolyData()
    vtk.vtkCommunicator.UnMarshalDataObject(dsa.numpyTovtkDataArray(buff), pd)
    print pd
'''

'''
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
   data = {'a': 7, 'b': 3.14}
   comm.send(data, dest=1, tag=11)
elif rank == 1:
   data = comm.recv(source=0, tag=11)
   print rank
   print data
'''   

'''
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()

if rank == 0:
   data = {'key1' : [7, 2.72, 2+3j],
           'key2' : ( 'abc', 'xyz')}
else:
   data = None
data = comm.bcast(data, root=0)

print rank
print data
'''


from mpi4py import MPI

comm = MPI.COMM_WORLD
size = comm.Get_size()
rank = comm.Get_rank()

if rank == 0:
   data = [(i+1)**2 for i in range(size)]
else:
   data = None
data = comm.scatter(data, root=0)

print 'r: ' + str(rank)
print data

assert data == (rank+1)**2




