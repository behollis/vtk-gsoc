import vtk
import h5py
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa
#import mayavi
#from mayavi import tools
#from mayavi import mlab
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors


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
    
def readStreamlines(x, y, f):
    slines = list()
    
    try:
        groups = f.keys()
        
        for mem in groups:
            dir = str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3)
        
            xlst = list(f[dir][0])
            ylst = list(f[dir][1])
            zlst = 0.0 * len(list(f[dir][0]))
        
            slines.append([xlst,ylst,zlst])
    except:
        print 'reading error for streamlines at: {0} {1}'.format(x,y)
        slines.append( [ [ERROR_CODE]*2,[ERROR_CODE]*2,[ERROR_CODE]*2 ] )
        
    #print slines
    return slines   
    
DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'
    
if __name__ == '__main__':
    
    f = h5py.File(DPATH+'0clusterFlowStruct.hdf5', 'r')
    f1 = h5py.File(DPATH+'1clusterFlowStruct.hdf5', 'r')
    f2 = h5py.File(DPATH+'2clusterFlowStruct.hdf5', 'r')
    
    
    sl_f1 = h5py.File(DPATH+'0.0.hdf5', 'r')
    sl_f2 = h5py.File(DPATH+'0.1.hdf5', 'r')
    sl_f3 = h5py.File(DPATH+'0.2.hdf5', 'r')
    sl_f4 = h5py.File(DPATH+'0.3.hdf5', 'r')
    sl_f5 = h5py.File(DPATH+'0.4.hdf5', 'r')
    sl_f6 = h5py.File(DPATH+'0.5.hdf5', 'r')
    
    slines1 = readStreamlines(35,24,sl_f1)
    #slines2 = readStreamlines(30,30,sl_f2)
    slines3 = readStreamlines(35,24,sl_f3)
    #slines4 = readStreamlines(30,30,sl_f4)
    slines5 = readStreamlines(35,24,sl_f5)
    #slines6 = readStreamlines(30,30,sl_f6)
    
    slines = slines1 + slines3 + slines5 

        
    
    s = np.zeros(shape=(125,125), dtype=float)
    
    for x in range(1,41):
        for y in range(0,125):
            s[y,x] = f['clusterCnt'][x,y] 
            #print 'x ' + str(x)
            #print 'y ' + str(y)
            #print s[y,x]
            #print f['clusterCnt'][x,y] 
    for x in range(41,81):
        for y in range(0,125):
            s[y,x] = f1['clusterCnt'][x,y] 
            #print 'x ' + str(x)
            #print 'y ' + str(y)
            #print s[y,x]
            #print f['clusterCnt'][x,y] 
            
    for x in range(82,125):
        for y in range(0,125):
            s[y,x-2] = f2['clusterCnt'][x,y] 
            #print 'x ' + str(x)
            #print 'y ' + str(y)
            #print s[y,x]
            #print f['clusterCnt'][x,y] 
            
    # make a color map of fixed colors
    cmap = colors.ListedColormap(['white', 'blue', 'green', 'red', 'orange', 'yellow'])
    bounds=[0,1,2,3,4,5]
    norm = colors.BoundaryNorm(bounds, cmap.N)
   
    imgplot = plt.imshow(s, interpolation='none', origin='upper',cmap=cmap, norm=norm)
    
    
    for sl in slines:
        plt.plot(sl[0], sl[1], color='black', linewidth=0.5, alpha=0.8)
        
       
    slines1 = readStreamlines(15,70,sl_f1)
    #slines2 = readStreamlines(30,30,sl_f2)
    slines3 = readStreamlines(15,70,sl_f3)
    #slines4 = readStreamlines(30,30,sl_f4)
    slines5 = readStreamlines(15,70,sl_f5)
    #slines6 = readStreamlines(30,30,sl_f6)
    
    slines = slines1 + slines3 + slines5 
    for sl in slines:
        plt.plot(sl[0], sl[1], color='black', linewidth=0.5, alpha=0.8)
    
        
    # make a color bar
    plt.colorbar(imgplot, cmap=cmap, norm=norm, boundaries=bounds, ticks=[-1,0,1,2,3,4,5])
    #imgplot.set_interpolation('none')
    #imgplot.set_cmap('hot')
    plt.gca().invert_yaxis()
    plt.show()
    
    
    
   
    