import vtk
import h5py
from vtk.util.vtkAlgorithm import VTKPythonAlgorithmBase
from vtk.numpy_interface import dataset_adapter as dsa
#import mayavi
#from mayavi import tools
#from mayavi import mlab
import numpy as np



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
    
    f = h5py.File(DPATH+'ftva.hdf5', 'r')
    
    '''
    s = np.ndarray(shape=(41,125,30))
    
    for x in range(0,41):
        for y in range(0,125):
            for m in range(0,30):
                s[x,y,m] = f['majEigenVal'][x,y,m] 
                #print f['majEigenVal'][x][y][m] 
                
    sf = mlab.pipeline.scalar_field(s)
    mlab.pipeline.volume(s)
    '''
        
    alg = HDF5Source()
    alg.SetFileName(DPATH + 'ftva.hdf5')
    alg.Update()
    
    
    im = alg.GetOutputPort()
    #m = vtk.vtkPolyDataMapper()
    #m.SetInputConnection(im)

    #cf = vtk.vtkContourFilter()
    #cf.SetInputConnection(alg.GetOutputPort())
    #cf.SetValue(0, 200)
    
    # The volume will be displayed by ray-cast alpha compositing.
    # A ray-cast mapper is needed to do the ray-casting, and a
    # compositing function is needed to do the compositing along the ray.
    #rayCastFunction = vtk.vtkVolumeRayCastCompositeFunction()
    
    #volumeMapper = vtk.vtkVolumeRayCastMapper()
    volumeMapper = vtk.vtkSmartVolumeMapper()
    volumeMapper.SetInputConnection( im )
    #volumeMapper.SetVolumeRayCastFunction(rayCastFunction)
    
    attr = f['majEigenVal'][:]
    min = attr.min()
    max = attr.max()
    
    # The color transfer function maps voxel intensities to colors.
    # It is modality-specific, and often anatomy-specific as well.
    # The goal is to one color for flesh (between 500 and 1000)
    # and another color for bone (1150 and over).
    volumeColor = vtk.vtkColorTransferFunction()
    volumeColor.AddRGBPoint(min,    0.0, 0.0, 0.0)
    #volumeColor.AddRGBPoint(500,  1.0, 0.5, 0.3)
    #volumeColor.AddRGBPoint(1000, 1.0, 0.5, 0.3)
    volumeColor.AddRGBPoint(max, 1.0, 0.0, 0.0)
    
    
    # The opacity transfer function is used to control the opacity
    # of different tissue types.
    volumeScalarOpacity = vtk.vtkPiecewiseFunction()
    volumeScalarOpacity.AddPoint(min,    0.00)
    #volumeScalarOpacity.AddPoint(500,  0.15)
    #volumeScalarOpacity.AddPoint(1000, 0.15)
    volumeScalarOpacity.AddPoint(max, 0.85)
    
    # The gradient opacity function is used to decrease the opacity
    # in the "flat" regions of the volume while maintaining the opacity
    # at the boundaries between tissue types.  The gradient is measured
    # as the amount by which the intensity changes over unit distance.
    # For most medical data, the unit distance is 1mm.
    #volumeGradientOpacity = vtk.vtkPiecewiseFunction()
    #volumeGradientOpacity.AddPoint(0,   0.0)
    #volumeGradientOpacity.AddPoint(90,  0.5)
    #volumeGradientOpacity.AddPoint(100, 1.0)
    
    # The VolumeProperty attaches the color and opacity functions to the
    # volume, and sets other volume properties.  The interpolation should
    # be set to linear to do a high-quality rendering.  The ShadeOn option
    # turns on directional lighting, which will usually enhance the
    # appearance of the volume and make it look more "3D".  However,
    # the quality of the shading depends on how accurately the gradient
    # of the volume can be calculated, and for noisy data the gradient
    # estimation will be very poor.  The impact of the shading can be
    # decreased by increasing the Ambient coefficient while decreasing
    # the Diffuse and Specular coefficient.  To increase the impact
    # of shading, decrease the Ambient and increase the Diffuse and Specular.
    volumeProperty = vtk.vtkVolumeProperty()
    volumeProperty.SetColor(volumeColor)
    volumeProperty.SetScalarOpacity(volumeScalarOpacity)
    #volumeProperty.SetGradientOpacity(volumeGradientOpacity)
    volumeProperty.SetInterpolationTypeToLinear()
    volumeProperty.ShadeOn()
    volumeProperty.SetAmbient(0.4)
    volumeProperty.SetDiffuse(0.6)
    volumeProperty.SetSpecular(0.2)

    
    # The vtkVolume is a vtkProp3D (like a vtkActor) and controls the position
    # and orientation of the volume in world coordinates.
    volume = vtk.vtkVolume()
    volume.SetMapper(volumeMapper)
    volume.SetProperty(volumeProperty)
    
    # With almost everything else ready, its time to initialize the renderer and window, as well as creating a method for exiting the application
    renderer = vtk.vtkRenderer()
    renderWin = vtk.vtkRenderWindow()
    renderWin.AddRenderer(renderer)
    renderInteractor = vtk.vtkRenderWindowInteractor()
    renderInteractor.SetRenderWindow(renderWin)
   
    # We add the volume to the renderer ...
    renderer.AddVolume(volume)
    # ... set background color to white ...
    renderer.SetBackground(1, 1, 1)
    # ... and set window size.
    renderWin.SetSize(400, 400)
   
    # A simple function to be called when the user decides to quit the application.
    def exitCheck(obj, event):
        if obj.GetEventPending() != 0:
            obj.SetAbortRender(1)
            
    # Tell the application to use the function as an exit check.
    renderWin.AddObserver("AbortCheckEvent", exitCheck)
   
    renderInteractor.Initialize()
    # Because nothing will be rendered without any input, we order the first render manually before control is handed over to the main-loop.
    
    vtkAxesActor_axes = vtk.vtkAxesActor()
 
    vtkOrientationMarkerWidget_widget = vtk.vtkOrientationMarkerWidget()
    vtkOrientationMarkerWidget_widget.SetOutlineColor( 0.9300, 0.5700, 0.1300 )
    vtkOrientationMarkerWidget_widget.SetOrientationMarker( vtkAxesActor_axes )
    vtkOrientationMarkerWidget_widget.SetInteractor( renderInteractor )
    vtkOrientationMarkerWidget_widget.SetViewport( 0.0, 0.0, 0.4, 0.4 )
    vtkOrientationMarkerWidget_widget.SetEnabled( 1 )
    vtkOrientationMarkerWidget_widget.InteractiveOn()
 
    renderer.ResetCamera()
    
    renderWin.Render()
    renderInteractor.Start()
    
    #while True:
    #    renWin.Render()
    
   
    