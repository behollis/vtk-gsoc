from sklearn.cluster import DBSCAN
import math
import h5py
import numpy as np

from mpi4py import MPI

DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'

import numpy.linalg as la

'''
implementation of entropy based streamline cluster from:

Chen, ChengKai, et al. An illustrative visualization framework for 3d vector fields. Computer Graphics Forum. 
Vol. 30. No. 7. Blackwell Publishing Ltd, 2011.
'''

def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def calcStreamlineFeatures(slines, entropy_only = False, sline_dict = None):
    slines_features = list()
    
    for sl in slines:
        
       if len(sl[0]) < 3 and len(sl[1]) < 3 and len(sl[2]) < 3:
           continue 
         
       slf = list() #list of streamline feature vectors to cluster
       #find first, middle, and last points on sl
       begin_pt = [ sl[0][0], sl[1][0] ] #, 0.0 ]
       half_num_pts = len(sl[0])/2
       
       # TODO: need to compute this based on arc length approx of sline
       mid_pt = [ sl[ 0 ][ half_num_pts ], sl[ 1 ][ half_num_pts ] ] #, 0.0  ]
       
       end_pt = [ sl[0][-1], sl[1][-1] ] #, 0.0 ]
       slf = begin_pt + mid_pt + end_pt 
       
       linear_entropy = 0.0
       angl_entropy = 0.0
       
       num_pts = len( sl[0] ) 
       num_segs = num_pts - 1
       
       #find total length of streamline and length of each seg
       sl_seg_lengths = list()
       total_len = 0.0
       for j in range( 0, num_pts - 1 ):
           pt0 = ( sl[0][j], sl[1][j], 0.0 )
           pt1 = ( sl[0][j+1], sl[1][j+1], 0.0 )
           
           dist = np.linalg.norm( np.array(pt1) - np.array(pt0) )
           sl_seg_lengths.append( dist )
           
           total_len += dist
           
       #find midpoint location of arc length
       mid_pt_idx = mid_pt
       arc_length = 0.0
       for j in range( 0, num_pts - 1 ):
           pt0 = ( sl[0][j], sl[1][j], 0.0 )
           pt1 = ( sl[0][j+1], sl[1][j+1], 0.0 )
           
           dist = np.linalg.norm( np.array(pt1) - np.array(pt0) )
           arc_length += dist
          
           if arc_length >= ( total_len / 2.0 ):
               mid_pt_idx = j
               break
           
       # reassign proper mid point    
       #print slf
       
       slf[2] = sl[ 0 ][ mid_pt_idx ]
       slf[3] = sl[ 1 ][ mid_pt_idx ]
       
       #slf[5] = 0.0  
       
       #print slf
       
       '''
       #summation of linear entropies
       sum_lentropy = 0.0
       for s in sl_seg_lengths:
           #sum_lentropy += s + ( math.log(s, 2) / total_len )
           sum_lentropy += ( s / total_len ) * ( math.log(s / total_len, 2) / total_len )
       
       #lentropy = -1.0 * ( sum_lentropy / ( math.log(total_len + 1,2) * total_len ) )
       lentropy = -1.0 * sum_lentropy / math.log( total_len + 1,2 )
       
       slf.append( lentropy )
       
       #find total angular variation and angle between each segment
       sl_seg_ang = list()
       total_ang = 0.0
       total_ang_change = 0.0
       for j in range( 1, num_pts - 1 ):
           pt0 = ( sl[0][j-1], sl[1][j-1], 0.0 )
           pt1 = ( sl[0][j]  , sl[1][j]  , 0.0 )
           pt2 = ( sl[0][j+1], sl[1][j+1], 0.0 )
           
           v0 = np.array(pt1) - np.array(pt0)
           v1 = np.array(pt2) - np.array(pt1)
           
           ang = py_ang(v0,v1)
           #print 'ang: ' + str(ang)
           
           sl_seg_ang.append(ang)
           
           total_ang += math.fabs(ang)
           total_ang_change += ang
           
       #summation of linear entropies
       sum_aentropy = 0.0
       for a in sl_seg_ang:
           if a != 0.0:
               #sum_aentropy += a + ( math.log( a, 2) / total_ang )
               sum_aentropy += ( a / total_ang ) * ( math.log(a / total_ang, 2) / total_ang )
               
           
       #aentropy = -1.0 * ( sum_aentropy / ( math.log(total_ang,2) * total_ang ) )
       aentropy = -1.0 * sum_aentropy / math.log( total_ang,2 )
       
       slf.append( aentropy )
       
       if entropy_only: 
           slines_features.append( [ slf[9], slf[10] ] )
       '''
        
       #else:
       slines_features.append( slf )
           
       #add entry to dictionary
       if sline_dict is not None:
           sline_dict[tuple(slf)] = sl
    
       #print 'l entropy: ' + str( slf[9] )
       #print 'a entropy: ' + str( slf[10] )
            
    # return feature vectors of streamlines
    return slines_features

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
   
NUM_PROC = 8 #number of processes / cores per node
NUM_NODES = 3
X_EXT = 125
Y_EXT = 125 
SLICE_Y = int( Y_EXT / float( NUM_PROC ) )
SLICE_X = int ( X_EXT / float( NUM_NODES ) )
STEPS = 500
INTERVAL = 10
MEMBERS = 1000
MINTERVAL = 10
ERROR_CODE = -1
        
if __name__ == '__main__':
    
    comm = MPI.COMM_WORLD
    rank = comm.rank
    
    #determine node
    node = 0#rank / NUM_PROC
    print 'node: ' + str(node)
    print 'rank: ' + str(rank)
    
    #print 'reading hdf5 file...'
    if node == 0: #node = 1
        f0 = h5py.File(DPATH+'0.0.hdf5', 'r', driver='mpio', comm=comm)
        f1 = h5py.File(DPATH+'0.1.hdf5', 'r',driver='mpio', comm=comm)
        f2 = h5py.File(DPATH+'0.2.hdf5', 'r',driver='mpio', comm=comm)
        f3 = h5py.File(DPATH+'0.3.hdf5', 'r',driver='mpio', comm=comm)
        f4 = h5py.File(DPATH+'0.4.hdf5', 'r',driver='mpio', comm=comm)
        f5 = h5py.File(DPATH+'0.5.hdf5', 'r',driver='mpio', comm=comm)
    elif node == 3: #node = 3
        f0 = h5py.File(DPATH+'2.12.hdf5', 'r', driver='mpio', comm=comm)
        f1 = h5py.File(DPATH+'2.13.hdf5', 'r',driver='mpio', comm=comm)
        f2 = h5py.File(DPATH+'2.14.hdf5', 'r',driver='mpio', comm=comm)
        f3 = h5py.File(DPATH+'2.15.hdf5', 'r',driver='mpio', comm=comm)
        f4 = h5py.File(DPATH+'2.16.hdf5', 'r',driver='mpio', comm=comm)
        f5 = h5py.File(DPATH+'2.17.hdf5', 'r',driver='mpio', comm=comm)
    else:
        f0 = h5py.File(DPATH+'1.0.hdf5', 'r', driver='mpio', comm=comm)
        f1 = h5py.File(DPATH+'1.1.hdf5', 'r',driver='mpio', comm=comm)
        f2 = h5py.File(DPATH+'1.2.hdf5', 'r',driver='mpio', comm=comm)
        f3 = h5py.File(DPATH+'1.3.hdf5', 'r',driver='mpio', comm=comm)
        f4 = h5py.File(DPATH+'1.4.hdf5', 'r',driver='mpio', comm=comm)
        f5 = h5py.File(DPATH+'1.5.hdf5', 'r',driver='mpio', comm=comm)
        
        
    '''
    f0.atomic = True
    f1.atomic = True
    f2.atomic = True
    f3.atomic = True
    f4.atomic = True
    f5.atomic = True
    '''   
        
    print 'finished reading hdf5 files'
        
    output = h5py.File( DPATH+str(node)+'clusterFlowStruct.hdf5', 'w', driver='mpio', comm=comm)
    
    dset = output.create_dataset( name='clusterCnt',shape = ( X_EXT, Y_EXT), dtype='f' )
                                  
    # streamline clustering for seed...
    print 'end x-range: {0}'.format((node+1)*SLICE_X) 
    for x in range( node*SLICE_X, (node+1)*SLICE_X ): 
        for y in range( rank*SLICE_Y, (rank+1)*SLICE_Y ):
            
            if x == 0:
                #avoid edge case
                dset[x,y] = ERROR_CODE  
                continue
            
            print 'x{0} y{1} on rank: {2}'.format(x,y,rank)
             
            slines0 = readStreamlines(x, y, f0)
            slines2 = readStreamlines(x, y, f2)
            slines4 = readStreamlines(x, y, f4)
            #slines1 = readStreamlines(x, y, f1)
            #slines3 = readStreamlines(x, y, f3)
            #slines5 = readStreamlines(x, y, f5)
            
            #sltotal = slines0 + slines1 + slines2 + slines3 + slines4 + slines5
            sltotal = slines0 + slines2 + slines4 
            
            
            #fig = plt.figure()
            #ax = fig.gca(projection='3d')
            #ax.scatter(29,30,0,c='black',s=30)
            #plt.title('All members, Size of Cluster: {0}'.format(len(sltotal)))
            #for l in sltotal:
            #    ax.plot(l[0], l[1], l[2])    
            #ax.legend()
            #plt.show()
            
            
            #calculate feature vectors for streamlines
            feat = calcStreamlineFeatures(sltotal)
            
            FEATURES = 2 * 3
            if len(feat) >= FEATURES:
            
                print 'getting features...'
                feat_total = np.array( feat )
                feat_total.reshape(len(feat_total), FEATURES )
                
                #dbscan on 
                db = DBSCAN(eps=5.0).fit(feat_total)
                
                core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
                core_samples_mask[db.core_sample_indices_] = True
                labels = db.labels_
                
                # Number of clusters in labels, ignoring noise if present.
                n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
                
                print 'xy on rank: ' + str(x) + ' ' + str(y) + ' ' + str(rank)
                print 'clusters: ' + str( n_clusters_ )
                
                # test to see if this is where hang is
                
                print 'writing to hdf5 file...'    
                dset[x,y] = n_clusters_  
                print 'finished writing to hdf5!'
                
            else:
                print 'xy on rank: ' + str(x) + ' ' + str(y) + ' ' + str(rank)
                print 'reading error'
                    
                dset[x,y] = ERROR_CODE
                
                print 'wrote to dset'
               
    #comm.Barrier() 
    print 'closing hdf5 files...'
    f0.close()
    print 'closed f0'
    f1.close()
    print 'closed f1'
    f2.close()
    print 'closed f2'
    f3.close()
    print 'closed f3'
    f4.close()
    print 'closed f4'
    f5.close()
    print 'closed f5'
    output.close()
    print 'output closed'
   
        
    
    print 'rank: ' + str(rank) + ' finished!'
    
    
    
    
   
        
    
    
   
    