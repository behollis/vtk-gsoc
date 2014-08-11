import numpy
import vtk
import os
import sys
from subprocess import *

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


def readStreamline(stfile,stlst_f, stlst_b):
      
    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(stfile)
    
    try:
        reader.Update()
    except:
        print 'Can\'t read vtp file.'
        return

    #one streamline per vtp
    try:
        vtkPolyData_vtp = reader.GetOutput()
    except:
        print 'Can\'t read vtp file.'
        return
    
    #forward sl
    #vtkPolyLine_sl = vtkPolyData_vtp.GetCell(0)
    #vtkPoints_pts = vtkPolyLine_sl.GetPoints()
    #stlst_f.append(vtkPoints_pts)
    
    #backward sl
    if vtkPolyData_vtp.GetNumberOfCells() >= 2:
        vtkPolyLine_sl = vtkPolyData_vtp.GetCell(1)
        vtkPoints_pts = vtkPolyLine_sl.GetPoints()
        
        stlst_b.append(vtkPoints_pts)
    else:
        pass
        #print 'missing backward streamline for: ' + str(stfile)

def collectEnsembleStreamlines(x, y, stlst_f, stlst_b):
    ''' 
    Read all points for each streamline thru a seed for the ensemble. 
    '''
    ls_output = Popen(['ls'], stdout=PIPE)
    dirs = ls_output.communicate()[0].split('\n')
    
    for mem_dir in dirs[0:MEMBERS]:
        d = PATH + mem_dir + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3) +'/' 
        
        try:
            mem = mem_dir[3:]
            path_file = d+'sline_M'+ mem +'_X'+str(x).zfill(3) + '_Y' + str(y).zfill(3) +'.vtp'
            #if os.path.isfile(path_file):
            readStreamline(path_file, stlst_f, stlst_b)
        except:
            print 'Error reading streamline: ' + str(mem_dir)
    
def performStats(x,y,st_f, st_b):
    print 'performing stats on: ' + str(x) + ', ' + str(y)
    
    dir = OUT_PATH + 'ftva/x' + str(x).zfill(3) + '/y' + str(y).zfill(3) +'/' 
    if not os.path.exists(dir):
        os.makedirs(dir)
  
    out = open(dir+'ftva.txt', 'w+')
    
    if len(st_b) > 2:
        for pt_id in range(0,STEPS,INTERVAL):
            #print 'sline point: ' + str(pt_id)
             
            comp1_array = vtk.vtkDoubleArray()
            comp1_array.SetNumberOfComponents(1)
            comp2_array = vtk.vtkDoubleArray()
            comp2_array.SetNumberOfComponents(1)
            
            for line in st_b:
                try:
                    pt_tuple = line.GetPoint(pt_id)
                    
                    comp1_array.SetName( 'x' )
                    comp1_array.InsertNextValue(pt_tuple[0])
                 
                    comp2_array.SetName( 'y' )
                    comp2_array.InsertNextValue(pt_tuple[1])
                except:
                    #print 'missing points on member streamline: ' + str(line)
                    continue
             
            try:       
                var = calcPCA(comp1_array, comp2_array)
            except:
                #print 'PCA couldn\'t be calculated.'
                var = -1 
   
            out.write(str(var)+'\n')
            
        out.close()
    
    else:
        #print 'PCA couldn\'t be calculated, due to insufficent streamlines.'
        var = -1 
    
    
def main():
    for x in range(X_STR+BEGIN_OFFSET, X_END-END_OFFSET, SEED_RES):
        for y in range(Y_STR+BEGIN_OFFSET, Y_END-END_OFFSET, SEED_RES):
            collectEnsembleStreamlines(x, y, SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            try:
                performStats(x,y,SEED_STREAMLINES_F, SEED_STREAMLINES_B)
            except:
                #pass
                print 'Error computing PCA for: ' + str(x) + ' ' + str(y)
            del SEED_STREAMLINES_F[:]; del SEED_STREAMLINES_B[:]
            
            
if __name__ == '__main__':
    PATH = '/home/data3/lockExSt/ts00050/'
    OUT_PATH = '/home/nfs2/lockExSt/ts00050/'
    SEED_STREAMLINES_F = list()
    SEED_STREAMLINES_B = list()
    
    STEPS = 1000 # number of points on streamlines
    MEMBERS = 999
    X_EXT = 25 * 5
    Y_EXT = 25 * 5
    NUM_CORES = 26 # need this to be a squared integer
    SEED_RES = 1 # regularity of seed sampling for FTVA
    INTERVAL = 1 # interval of integration steps for multi-step FTVA
    
    wc = vtk.vtkMPIController()
    gsize = wc.GetNumberOfProcesses()
    grank = wc.GetLocalProcessId()
    if gsize % 2 != 0:
        if grank == 0:
            print 'Only an even number of ranks is support'
        sys.exit(0)
    
    localSize = gsize / NUM_CORES
    localGroup = grank / localSize
    
    #DEBUG, run as single process
    #grank = 0
    
    if grank == 25:
        print 'grank exit for: ' + str(grank)
        exit()
        
    os.chdir(PATH)
    
    # initialize offsets
    BEGIN_OFFSET = 0
    END_OFFSET = 0
    X_STR = 0; Y_STR = 0
    X_END = X_EXT; Y_END = Y_EXT 
    
    DIV = 5#int(NUM_CORES**0.5) # index for field block along each dimension
    SLICE = X_EXT / DIV # x or y length in field dimension
    
    mod_x = grank % DIV # determine index along x dir for block
    div_y = grank / DIV # determine index along y dir for block
    
    X_STR = mod_x * SLICE
    X_END = (mod_x + 1) * SLICE
    Y_STR = div_y * SLICE
    Y_END = (div_y + 1) * SLICE
    
    print 'grank starting: ' + str(grank)
    print '\tX_STR: ' + str(X_STR)
    print '\tX_END: ' + str(X_END)
    print '\tY_STR: ' + str(Y_STR)
    print '\tY_END: ' + str(Y_END)
        
    main()
    
    print 'grank finished: ' + str(grank)
    