from mpi4py import MPI
import h5py

DPATH = '/home/behollis/DATA/out/lockExSt/ts00050/'


f0 = h5py.File(DPATH+'0.0.hdf5', 'r', driver='mpio', comm=MPI.COMM_WORLD)
rank = MPI.COMM_WORLD.rank  # The process ID (integer 0-3 for 4-process run)

f = h5py.File('parallel_test22.hdf5', 'w', driver='mpio', comm=MPI.COMM_WORLD)

dset = f.create_dataset('test', (4,), dtype='i')
dset[rank] = rank + 8

f.close()