import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
import math
import random


DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'

'''
implementation of entropy based streamline cluster from:

Chen, ChengKai, et al. An illustrative visualization framework for 3d vector fields. Computer Graphics Forum. 
Vol. 30. No. 7. Blackwell Publishing Ltd, 2011.
'''

import numpy.linalg as la
 
def py_ang(v1, v2):
    """ Returns the angle in radians between vectors 'v1' and 'v2'    """
    cosang = np.dot(v1, v2)
    sinang = la.norm(np.cross(v1, v2))
    return np.arctan2(sinang, cosang)

def calcStreamlineFeatures(slines):
    slines_features = list()
    
    for sl in slines:
        
       if len(sl[0]) < 3 and len(sl[1]) < 3 and len(sl[2]) < 3:
           continue 
         
       slf = list() #list of streamline feature vectors to cluster
       #find first, middle, and last points on sl
       begin_pt = [ sl[0][0], sl[1][0], 0.0 ]
       half_num_pts = len(sl[0])/2
       
       # TODO: need to compute this based on arc length approx of sline
       mid_pt = [ sl[ 0 ][ half_num_pts ], sl[ 1 ][ half_num_pts ], 0.0  ]
       
       end_pt = [ sl[0][-1], sl[1][-1], 0.0 ]
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
       slf[3] = sl[ 0 ][ mid_pt_idx ]
       slf[4] = sl[ 1 ][ mid_pt_idx ]
       slf[5] = 0.0  
       
       #print slf
       
       
       #summation of linear entropies
       sum_lentropy = 0.0
       for s in sl_seg_lengths:
           sum_lentropy += s + ( math.log(s, 2) / total_len )
       
       lentropy = -1.0 * ( sum_lentropy / ( math.log(total_len + 1,2) * total_len ) )
       
       slf.append( lentropy )
       
       #find total angular variation and angle between each segment
       sl_seg_ang = list()
       total_ang = 0.0
       for j in range( 1, num_pts - 1 ):
           pt0 = ( sl[0][j-1], sl[1][j-1], 0.0 )
           pt1 = ( sl[0][j]  , sl[1][j]  , 0.0 )
           pt2 = ( sl[0][j+1], sl[1][j+1], 0.0 )
           
           v0 = np.array(pt1) - np.array(pt0)
           v1 = np.array(pt2) - np.array(pt1)
           
           ang = py_ang(v0,v1)
           sl_seg_ang.append(ang)
           
           total_ang += math.fabs(ang)
           
       #summation of linear entropies
       sum_aentropy = 0.0
       for a in sl_seg_ang:
           if a != 0.0:
               sum_aentropy += a + ( math.log( a, 2) / total_ang )
           
       aentropy = -1.0 * ( sum_aentropy / ( math.log(total_ang,2) * total_ang ) )
       
       slf.append( aentropy )
       
        
       slines_features.append( slf )
       
    # return feature vectors of streamlines
    return slines_features
        
def readStreamlines(x, y, f, mem_min=None, mem_max=None, ext_xmax = None, \
                    ext_ymax = None, ext_xmin = None, ext_ymin = None):
    ''' loads streamlines thru seed with spatial extents / member extents '''
    
    slines = list()
    xlst =list()
    ylst = list()
    zlst = list()
    groups = f.keys()
    
    for mem in groups:
        dir = str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3)
    
        mem_n = int(mem[3:]) 
        if (mem_min == None and mem_max == None) or (mem_n >= mem_min and mem_n < mem_max):
            try:
                xlst = list(f[dir][0])
                ylst = list(f[dir][1])
                zlst = 0.0 * len(list(f[dir][0]))
            except:
                print 'trying to read MEMBER or SEED that is not in this hdf5 file...'
        else:
            continue
            
        if len(xlst) > 0 and len(ylst) > 0:
            #collecting streamlines to cluster in region
            region_xlst = list()
            region_ylst = list()
            region_zlst = list()
            
            for idx in range(0, len(xlst)):
                x = xlst[idx]
                y = ylst[idx]
                if ( x <= ext_xmax ) and ( x >= ext_xmin ) \
                    and (y <= ext_ymax) and (y >= ext_ymin):
                    region_xlst.append(x)
                    region_ylst.append(y)
                    region_zlst.append(0.0)
                            
            if len(region_xlst) > 0:          
                slines.append([region_xlst, region_ylst, region_zlst])
                            
                    
    return slines

def clusterCrispFieldSlines(vec_field):
    #calculate feature vectors for streamlines
    feat = list()
    for sl in vec_field:
        feat_sl = calcStreamlineFeatures( sl )
        if len( feat_sl ) > 0:
            feat.append( feat_sl[ 0 ] )
    
    feat_total = np.array( feat )
    
    feat_total.reshape(len(feat_total), 3*3 + 2 )
    
    #dbscan on 
    db = DBSCAN(eps=50.0).fit(feat_total)
    
    return db
        
if __name__ == '__main__':
    #print 'reading hdf5 file...'
    f0 = h5py.File(DPATH+'0.0.hdf5', 'r')
    f1 = h5py.File(DPATH+'0.1.hdf5', 'r')
    f2 = h5py.File(DPATH+'0.2.hdf5', 'r')
    f3 = h5py.File(DPATH+'0.3.hdf5', 'r')
    f4 = h5py.File(DPATH+'0.4.hdf5', 'r')
    f5 = h5py.File(DPATH+'0.5.hdf5', 'r')
    
    
    print 'finished reading h5py file!'
    
    files = (f0,f1,f2,f3,f4,f5)
    
    ext_xmin = 5; ext_xmax = 10; ext_ymin = 5; ext_ymax = 15
    mem_min = 1; mem_max = 5
    # readstreamlines for region
    
    member_slines = list()
    SKIP = 1
    for member in range(mem_min, mem_max + 1, 1):
        region_slines = list()
        print member
        for f in files:
            for seed_y in range(ext_ymin, ext_ymax, SKIP):
                for seed_x in range(ext_xmin, ext_xmax, SKIP): 
                    #print 'reading {0} {1}'.format(seed_x, seed_y)
                    region_slines.append( readStreamlines(seed_x, seed_y, f, member, member+1, \
                                                         ext_xmax, ext_ymax, ext_xmin, ext_ymin) )
        member_slines.append( region_slines )
            
    
    for idx in range(len(member_slines)):
        fig = plt.figure()
        #ax = fig.gca(projection='3d')
        
        random.seed()
        c = (random.uniform(0.,1.), random.uniform(0.,1.), random.uniform(0.,1.))
        for g in member_slines[idx]:
            if len(g) != 0:
                for l in g:
                    #ax.plot(l[0], l[1], 0.0, color = c, linewidth = 0.3)
                    plt.plot(l[0], l[1], color = c, linewidth = 0.3)
                    #print 'color{0}'.format(c) 
                    #print l[0]
                    #print l[1]
                    #print l[2] 
        plt.title('Region Streamlines, member {0}'.format(idx))  
        #ax.legend()
        plt.savefig(DPATH+'regionslines{0}'.format(idx))
        
    for idx, m in enumerate(member_slines):
        db = clusterCrispFieldSlines(m)
        
        for c in set(db.labels_):
            fig = plt.figure()
            #ax = fig.gca(projection='3d')
            
            size_cluster = 0
            for l in range(0, len( db.labels_ )):
                if c == int(db.labels_[l] ):
                    if m[l][0][0] > 2 and m[l][0][1] > 2 and m[l][0][2] > 2:
                        size_cluster += 1
                        plt.plot( m[l][0][0], m[l][0][1] ) 
                    
            plt.title('Member: {0}, Label: {1}, Size: {2}'.format(idx,c,size_cluster))   
            #ax.legend()
            plt.show()
    
    '''
    # streamline clustering for seed...
    
    slines0 = readStreamlines(29,30, f0)
    slines2 = readStreamlines(29,30, f2)
    slines4 = readStreamlines(29,30, f4)
    slines1 = readStreamlines(29,30, f1)
    slines3 = readStreamlines(29,30, f3)
    slines5 = readStreamlines(29,30, f5)
    
    sltotal = slines0 + slines1 + slines2 + slines3 + slines4 + slines5
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.scatter(29,30,0,c='black',s=30)
    plt.title('All members, Size of Cluster: {0}'.format(len(sltotal)))
    for l in sltotal:
        ax.plot(l[0], l[1], l[2])    
    
    ax.legend()
    plt.show()
    
    
    #calculate feature vectors for streamlines
    feat = calcStreamlineFeatures(sltotal)
    #feat2 = calcStreamlineFeatures(slines2)
    #feat4 = calcStreamlineFeatures(slines4)
    
    feat_total = np.array( feat )#feat0 + feat2 + feat4 )
    feat_total.reshape(len(feat_total), 3*3 + 2 )
    
    #dbscan on 
    db = DBSCAN(eps=5.0).fit(feat_total)
    
    
    
    core_samples_mask = np.zeros_like(db.labels_, dtype=bool)
    core_samples_mask[db.core_sample_indices_] = True
    labels = db.labels_
    
    # Number of clusters in labels, ignoring noise if present.
    n_clusters_ = len(set(labels)) - (1 if -1 in labels else 0)
    
    
    for c in set(labels):
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        ax.scatter(29,30,0,c='black',s=30)
        
        size_cluster = 0
        for l in range(0, len(sltotal)):
            if c == int(labels[l]):
                size_cluster += 1
                ax.plot(sltotal[l][0], sltotal[l][1], sltotal[l][2]) 
                
        plt.title('Label: {0}, Size: {1}'.format(c,size_cluster))   
        
        ax.legend()
        plt.show()
    '''
    
    f0.close()
    f1.close() 
    f2.close()
    f3.close()   
    f4.close()
    f5.close()
    
    print 'finished!'
    
    
    
    
   
    