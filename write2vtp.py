import string
import numpy
import vtk

#http://www.vtk.org/Wiki/VTK/Writing_VTK_files_using_python
#http://www.shocksolution.com/microfluidics-and-biotechnology/visualization/python-vtk-paraview/
#http://docs.enthought.com/mayavi/mayavi/data.html#vtk-data-structures

EXT_X = 10
EXT_Y = 10
MEMBERS = 20
DATAPOINTS = EXT_X * EXT_Y

ftva = numpy.ndarray(shape=(EXT_X, EXT_Y))
i = vtk.vtkImageData()
i.SetDimensions(EXT_X,EXT_Y,1)
i.AllocateScalars(vtk.VTK_DOUBLE, 1)

for x in range(0,EXT_X):
    for y in range(0,EXT_Y):
        f = open('./out/'+str(x)+'_'+str(y)+'_eigenvalue.txt', 'r')
        ftv = float(string.strip(f.readline()))
        i.SetScalarComponentFromDouble(x, y, 0, 0, ftv)

w = vtk.vtkXMLImageDataWriter()    
w.SetInputData(i)             
w.SetFileName('ftvaField.vtp')         
w.Write()


    
    
    
    
    