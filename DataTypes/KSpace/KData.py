import numpy as np
import h5py
from scipy.io import loadmat, savemat
from DataTypes.Geometry.Patterns import InterleavedMatLabPattern

class KData(InterleavedMatLabPattern):
    """Pattern that also holds the data for each points"""
    def __init__(self, kpath, kdatapath, time_per_acquisition=None, alternated_spokes=True):
        super().__init__(path=kpath, time_per_acquisition=time_per_acquisition, alternated_spokes=alternated_spokes)
        with h5py.File(kdatapath,'r') as f:
            data = f['kdata'][:]
        shape = data.shape
        n_coils = shape[0]
        self.n_coils = n_coils
        n_readouts_per_spoke = shape[3]
        n_spokes_per_interleaf = shape[2]
        n_interleaves = shape[1]
        if not (n_readouts_per_spoke,n_spokes_per_interleaf,n_interleaves) == (self.n_readouts_per_spoke,self.n_spokes_per_interleaf,self.n_interleaves):
            raise ValueError("k-space file does not have the same characteristics as kdata file.")
        self.data = np.zeros((n_readouts_per_spoke,n_spokes_per_interleaf,n_interleaves,n_coils),dtype=np.complex)
        for i_readout in range(n_readouts_per_spoke):
            for i_spoke in range(n_spokes_per_interleaf):
                for i_interleaf in range(n_interleaves):
                    self.points[i_readout,i_spoke,i_interleaf].data = {}
                    for i_coil in range(n_coils):
                        self.data[i_readout,i_spoke,i_interleaf,i_coil] = np.complex(data[i_coil][i_interleaf][i_spoke][i_readout][0], data[i_coil][i_interleaf][i_spoke][i_readout][1])
                        self.points[i_readout,i_spoke,i_interleaf].data[i_coil] = self.data[i_readout,i_spoke,i_interleaf,i_coil]
        del(data)

    def sorting_indexes(self, func):
        sorted_points = np.vectorize(func)(self.points)
        indexes = np.dstack(np.unravel_index(np.argsort(sorted_points.ravel()), sorted_points.shape))
        return indexes[0]

    def time_sorted(self):
        indexes = self.sorting_indexes(lambda point: point.t)
        points = np.zeros((self.n_coils,self.n_points))
        for i_coil in range(self.n_coils):
            for index in indexes:
                ind = index[0]*self.n_readouts_per_spoke*self.n_spokes_per_interleaf + index[1]*self.n_readouts_per_spoke + index[2]
                points[i_coil][ind] = float(np.absolute(self.data[:,:,:,i_coil][tuple(index)]))
        return points
