import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from sklearn.cluster import DBSCAN
import math


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
        
def readStreamlines(x, y, f):
    slines = list()
    groups = f.keys()
    
    for mem in groups:
        dir = str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3)
    
        xlst = list(f[dir][0])
        ylst = list(f[dir][1])
        zlst = 0.0 * len(list(f[dir][0]))
    
        slines.append([xlst,ylst,zlst])
        
    return slines
        
if __name__ == '__main__':
    #print 'reading hdf5 file...'
    f0 = h5py.File(DPATH+'0.0.hdf5', 'r')
    f1 = h5py.File(DPATH+'0.1.hdf5', 'r')
    f2 = h5py.File(DPATH+'0.2.hdf5', 'r')
    f3 = h5py.File(DPATH+'0.3.hdf5', 'r')
    f4 = h5py.File(DPATH+'0.4.hdf5', 'r')
    f5 = h5py.File(DPATH+'0.5.hdf5', 'r')
    #print 'finished reading h5py file!'
    
   
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
    

    #print('Estimated number of clusters: %d' % n_clusters_)
    '''
    print("Homogeneity: %0.3f" % metrics.homogeneity_score(labels_true, labels))
    print("Completeness: %0.3f" % metrics.completeness_score(labels_true, labels))
    print("V-measure: %0.3f" % metrics.v_measure_score(labels_true, labels))
    print("Adjusted Rand Index: %0.3f"
          % metrics.adjusted_rand_score(labels_true, labels))
    print("Adjusted Mutual Information: %0.3f"
          % metrics.adjusted_mutual_info_score(labels_true, labels))
    print("Silhouette Coefficient: %0.3f"
          % metrics.silhouette_score(X, labels))
    '''
    
    
    '''
    ##############################################################################
    # Plot result
    import matplotlib.pyplot as plt
    
    # Black removed and is used for noise instead.
    unique_labels = set(labels)
    colors = plt.cm.Spectral(np.linspace(0, 1, len(unique_labels)))
    for k, col in zip(unique_labels, colors):
        if k == -1:
            # Black used for noise.
            col = 'k'
    
        class_member_mask = (labels == k)
    
        xy = feat_total[class_member_mask & core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=14)
    
        xy = feat_total[class_member_mask & ~core_samples_mask]
        plt.plot(xy[:, 0], xy[:, 1], 'o', markerfacecolor=col,
                 markeredgecolor='k', markersize=6)
    
    plt.title('Estimated number of clusters: %d' % n_clusters_)
    plt.show()
    '''
       
   
    f0.close(); f2.close(); f4.close()
    
    print 'finished!'
    
    
    
    
   
    