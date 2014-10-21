import h5py
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


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
       slf = [ sl[0], sl[-1], sl[ len(sl)/2 ] ]
       
       linear_entropy = 0.0
       angl_entropy = 0.0
       
       num_pts = len( sl ) 
       num_segs = num_pts - 1
       
       #find total length of streamline and length of each seg
       sl_seg_lengths = list()
       total_len = 0.0
       for j in range( 0, num_pts ):
           dist = np.linalg.norm(np.ndarray(sl[j])-np.ndarray(sl[j+1]))
           sl_seg_lengths.append(dist)
           total_len += dist
       
       #summation of linear entropies
       sum_lentropy = 0.0
       for s in sl_seg_lengths:
           sum_lentropy += s + ( math.log(s, 2) / total_len )
       
       lentropy = -1.0 * ( sum_lentropy / ( math.log(total_len + 1,2) * total_len ) )
       
       slf.append( lentropy )
       
       #find total angular variation and angle between each segment
       sl_seg_ang = list()
       total_ang = 0.0
       for j in range( 1, num_pts ):
           v0 = np.ndarray(sl[j]) - np.ndarray(sl[j-1])
           v1 = np.ndarray(sl[j+1]) - np.ndarray(sl[j])
           ang = py_ang(v0,v1)
           sl_seg_ang.append(ang)
           total_ang += math.abs(ang)
           
       #summation of linear entropies
       sum_aentropy = 0.0
       for a in sl_seg_ang:
           sum_aentropy += a + ( math.log(a, 2) / total_ang )
       
       aentropy = -1.0 * ( sum_aentropy / ( math.log(total_ang,2) * total_ang ) )
       
       slf.append( aentropy )
       
       sline_features.append( slf )
       
    # return feature vectors of streamlines
    return sline_features
        
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
    #f1 = h5py.File(DPATH+'0.1.hdf5', 'r')
    f2 = h5py.File(DPATH+'0.2.hdf5', 'r')
    #f3 = h5py.File(DPATH+'0.3.hdf5', 'r')
    f4 = h5py.File(DPATH+'0.4.hdf5', 'r')
    #f5 = h5py.File(DPATH+'0.5.hdf5', 'r')
    #print 'finished reading h5py file!'
    
   
    slines0 = readStreamlines(29,30, f0)
    slines2 = readStreamlines(29,30, f2)
    slines4 = readStreamlines(29,30, f4)
    #slines1 = readStreamlines(28,29, f1)
    #slines3 = readStreamlines(28,29, f3)
    #slines5 = readStreamlines(28,29, f5)
    
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    for l in slines0:
        ax.plot(l[0], l[1], l[2])
    #for l in slines1:
    #    ax.plot(l[0], l[1], l[2])
    for l in slines2:
        ax.plot(l[0], l[1], l[2])
    #for l in slines3:
    #    ax.plot(l[0], l[1], l[2])
    for l in slines4:
        ax.plot(l[0], l[1], l[2])
    #for l in slines5:
    #    ax.plot(l[0], l[1], l[2])
        
        
        
    ax.legend()
    plt.show()
   
   
    f0.close(); f2.close(); f4.close()
    
    print 'finished!'
    
    
    
    
   
    