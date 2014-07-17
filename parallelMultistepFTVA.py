import numpy
import vtk

def calcPCA(xarray, yarray):
    ''' Returns the eigenvalue of the covariance matrix for the terminal particle 
        positions. 
    '''
    
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
     
    eigenvalues = vtk.vtkDoubleArray()
    pcaStatistics.GetEigenvalues(eigenvalues)
    
    return eigenvalues.GetValue(0)

if __name__ == '__main__':
    
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName('/home/behollis/lockExchangeStreamlinesTs0050/0/1Bslines0000.vtp');
    reader.Update();
    
    vtkPolyData_vtp = reader.GetOutput()
    sl_count = vtkPolyData_vtp.GetNumberOfCells()
    
    for sl_id in range(0,sl_count):
        vtkPolyLine_sl = vtkPolyData_vtp.GetCell(sl_id)
        vtkPoints_sl = vtkPolyLine_sl.GetPoints()
        
        # verify seed 
        print vtkPoints_sl.GetPoint(0)
    
    vtkPolyLine_sl = vtkPolyData_vtp.GetCell(10)
    vtkPoints_pts = vtkPolyLine_sl.GetPoints()
    sl_num_pts = vtkPoints_pts.GetNumberOfPoints()
    
    for pt_id in range(0,sl_num_pts):
        pt_tuple = vtkPoints_pts.GetPoint(pt_id)
        print 'streamline: ' + str(10) + ' point id: ' + str(pt_id) + ': position: ' + str(pt_tuple)
    
    
    
       
