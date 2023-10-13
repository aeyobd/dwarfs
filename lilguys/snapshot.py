import h5py
import numpy as np


class Snapshot:
    _file = None

    def __init__(self, positions, velocities, 
            ids=None, potential=None,
            header = {}):
        self.pos = positions
        self.vel = velocities
        self.ids = ids
        self.potential = potential
        self.header = header

    @classmethod
    def file(cls, filename):
        f = h5py.File(filename, "r")
        pos = get_h5_vector(f, "Coordinates")
        vel = get_h5_vector(f, "Velocities")
        IDs = get_h5_vector(f, "ParticleIDs")
        acc = get_h5_vector(f, "Acceleration")
        header = get_header(f)
        return cls(pos, vel, IDs,
                header=header)


    def save(self, filename=None):
        pass

    def __del__(self):
        if self._file:
            self._file.close()

    def __str__(self):
        pass

    def __len__(self):
        pass

    @property
    def x(self):
        return self.pos[:,0]

    @property
    def y(self):
        return self.pos[:,1]

    @property
    def z(self):
        return self.pos[:,2]

    @property
    def m(self):
        return self.header["MassTable"][1]

def get_header(h5_f):
    h5_header = h5_f["Header"].attrs
    header = {}
    for key, val in h5_header.items():
        header[key] = val

    return header

def get_h5_vector(h5_f, key):
    return np.array(h5_f["PartType1/" + key][()])

def create_h5_vector(h5_f, key, val):
    h5_f.create_dataset("PartType1/" + key, data=val)


def set_h5_header(h5_f, key):
    pass
