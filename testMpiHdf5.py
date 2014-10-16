import numpy as np
import h5py
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.rank

DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'

print 'reading hdf5 file...'
f0 = h5py.File(DPATH+'0.0.hdf5', driver='mpio', comm=comm)
#f1 = h5py.File(DPATH+'0.1.hdf5', driver='mpio', comm=comm)
#f2 = h5py.File(DPATH+'0.2.hdf5', driver='mpio', comm=comm)
#f3 = h5py.File(DPATH+'0.3.hdf5', driver='mpio', comm=comm)
#f4 = h5py.File(DPATH+'0.4.hdf5', driver='mpio', comm=comm)
#f5 = h5py.File(DPATH+'0.5.hdf5', driver='mpio', comm=comm)
print 'finished reading h5py file!'

NUM_PROC = 6
Y_EXT = int ( 125.0 / float(NUM_PROC) )
SLICE = Y_EXT / NUM_PROC
X_EXT = 40

'''
# print seed and next location
for x in range( 0, X_EXT):
    for y in range( rank*SLICE +1, (rank+1)*SLICE ):
        for mem in range( 5, 150):
            dir = '/mem' + str(mem).zfill(4) + '/x' + str(x).zfill(3) + '/y' + str(y).zfill(3)
            
            pt0_0 = f0[dir][0][0]
            pt1_0 = f0[dir][1][0]
           
            
            print 'Rank: ' + str(rank)
            print 'first  streamline point: ' + str(pt0_0)
            print 'second streamline point: ' + str(pt1_0)
             
print 'Rank: ' + str(rank) + ' completed!'
'''
    