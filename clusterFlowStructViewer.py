import vtk
import h5py
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa
#import mayavi
#from mayavi import tools
#from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt


class MyInteractorStyle(vtk.vtkInteractorStyleTrackballCamera):
 
    def __init__(self,parent=None):
        self.AddObserver("MiddleButtonPressEvent",self.middleButtonPressEvent)
        self.AddObserver("MiddleButtonReleaseEvent",self.middleButtonReleaseEvent)
 
    def middleButtonPressEvent(self,obj,event):
        print "Middle Button pressed"
        self.OnMiddleButtonDown()
        return
 
    def middleButtonReleaseEvent(self,obj,event):
        print "Middle Button released"
        self.OnMiddleButtonUp()
        return
 

class HDF5Source(VTKPythonAlgorithmBase):
    def __init__(self):
        VTKPythonAlgorithmBase.__init__(self,
            nInputPorts=0,
            nOutputPorts=1, outputType='vtkImageData')

        self.__FileName = ""

    def RequestData(self, request, inInfo, outInfo):
        f = h5py.File(self.__FileName, 'r')
        data = f['majEigenVal'][:]
        output = dsa.WrapDataObject(vtk.vtkImageData.GetData(outInfo))
        # Note that we flip the dimensions here because
        # VTK's order is Fortran whereas h5py writes in
        # C order.
        output.SetDimensions(data.shape[::-1])
        output.PointData.append(data.ravel(), 'majEigenVal')
        output.PointData.SetActiveScalars('majEigenVal')
        return 1
    
    def RequestInformation(self, request, inInfo, outInfo):
        f = h5py.File(self.__FileName, 'r')
        # Note that we flip the shape because VTK is Fortran order
        # whereas h5py reads in C order. When writing we pretend that the
        # data was C order so we have to flip the extents/dimensions.
        dims = f['majEigenVal'].shape[::-1]
        info = outInfo.GetInformationObject(0)
        info.Set(vtk.vtkStreamingDemandDrivenPipeline.WHOLE_EXTENT(),
            (0, dims[0]-1, 0, dims[1]-1, 0, dims[2]-1), 6)
       
        return 1

    def SetFileName(self, fname):
        if fname != self.__FileName:
            self.Modified()
            self.__FileName = fname

    def GetFileName(self):
        return self.__FileName    
    
DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'
    
if __name__ == '__main__':
    
    f = h5py.File(DPATH+'0clusterFlowStruct.hdf5', 'r')
    
    s = np.ndarray(shape=(120,10))
    
    for x in range(1,10):
        for y in range(0,119):
            s[y,x] = f['clusterCnt'][x,y] 
            #print 'x ' + str(x)
            #print 'y ' + str(y)
            #print s[y,x]
            #print f['clusterCnt'][x,y] 
    
    imgplot = plt.imshow(s)
    plt.colorbar()
    #imgplot.set_interpolation('none')
    imgplot.set_cmap('hot')
    plt.gca().invert_yaxis()
    plt.show()
    
    
    
   
    