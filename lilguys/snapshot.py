import h5py
import numpy as np
from os import path


class Snapshot:
    _filename = None
    gadget4 = False

    def __init__(self, positions, velocities, 
            IDs=None, potential=None,
            header = {}):
        self.pos = positions
        self.vel = velocities
        if IDs is None:
            IDs = np.arange(len(self))
        self.IDs = IDs
        if potential is None:
            potential = np.zeros(len(self))
        self.potential = potential
        self.header = header

    @classmethod
    def file(cls, filename, gadget4=False):
        if gadget4:
            with h5py.File(filename, "r") as f:
                pos = get_h5_vector(f, "Position")
                vel = get_h5_vector(f, "Velocity")
                IDs = get_h5_vector(f, "ParticleIDs")
                header = get_h5_header(f)
            c = cls(pos, vel, IDs, header=header)
            c._filename = filename
            c.gadget4 = True
            return c

        else:
            with h5py.File(filename, "r") as f:
                pos = get_h5_vector(f, "Coordinates")
                vel = get_h5_vector(f, "Velocities")
                IDs = get_h5_vector(f, "ParticleIDs")
                acc = get_h5_vector(f, "Acceleration")
                pot = get_h5_vector(f, "Potential")
                header = get_h5_header(f)

        c = cls(pos, vel, IDs, header=header, potential=pot)
        c._filename = filename
        return c


    def save(self, filename=None):
        if filename is None and path.exists(self._filename):
            self._save()
            print("snapshot saved")
        else:
            if filename is not None:
                self._filename = filename
            self._save_new()
            print("snapshot saved at ", self._filename)
        return True

    def filter(self, filt, inplace=False):
        assert len(filt) == len(self)
        copy = self.copy()
        copy.pos = copy.pos[filt]
        copy.vel = copy.vel[filt]
        copy.IDs = copy.IDs[filt]
        if copy.potential is not None:
            copy.potential = copy.potential[filt]
        return copy

    def _save(self):

        with h5py.File(self._filename, "a") as f:
            if self.gadget4:
                set_h5_vector(f, "Position", self.pos)
                set_h5_vector(f, "Velocity", self.vel)
            else:
                set_h5_vector(f, "Coordinates", self.pos)
                set_h5_vector(f, "Velocities", self.vel)

            set_h5_vector(f, "ParticleIDs", self.IDs)
            set_h5_header(f, self.header)

    def _save_new(self):
        with h5py.File(self._filename, "w") as f:
            part = f.create_group("PartType1")
            f.create_group("Header")

            set_h5_vector(f, "ParticleIDs", self.IDs)
            if self.gadget4:
                set_h5_vector(f, "Position", self.pos)
                set_h5_vector(f, "Velocity", self.vel)
            else:
                set_h5_vector(f, "Coordinates", self.pos)
                set_h5_vector(f, "Velocities", self.vel)
                set_h5_vector(f, "Potential", self.potential)
                set_h5_vector(f, "PertAccel", np.zeros(len(self)))
                set_h5_vector(f, "Acceleration", np.zeros(len(self)))
                set_h5_vector(f, "TimeStep", np.zeros(len(self)))
            # set_h5_vector(f, "Acceleration", self.acc)
            set_h5_header(f, self.header)


    def shift(self, p0, v0=[0,0,0], inplace=False):
        if inplace:
            self.pos += p0
            self.vel += v0
            return self
        else:
            shifted = self.copy()
            shifted.pos += p0
            shifted.vel += v0
            return shifted

    def scale(self, r_scale=1, v_scale=1, m_scale=1, inplace=False):
        if inplace:
            self.pos *= r_scale
            self.vel *= v_scale
            self.m *= m_scale
        else:
            scaled = self.copy()
            scaled.pos *= r_scale
            scaled.vel *= v_scale
            scaled.m *= m_scale
            return scaled


    def copy(self):
        pos = np.copy(self.pos)
        vel = np.copy(self.vel)
        IDs = np.copy(self.IDs)
        potential = np.copy(self.potential)
        header = self.header.copy()
        return Snapshot(pos, vel, IDs=IDs, potential=potential, header=header)


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
    def v_x(self):
        return self.vel[:,0]

    @property
    def v_y(self):
        return self.vel[:,1]

    @property
    def v_z(self):
        return self.vel[:,2]

    @property
    def r(self):
        return np.sqrt(self.x**2 + self.y**2 + self.z**2)

    @property
    def v(self):
        return np.sqrt(self.v_x**2 + self.v_y**2 + self.v_z**2)


    @property
    def m(self):
        return self.header["MassTable"][1]

    @m.setter
    def m(self, a):
        self.header["MassTable"][1] = a

    def __str__(self):
        return f"snapshot with {len(self)} particles"

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.x)




def get_h5_header(h5_f):
    h5_header = h5_f["Header"].attrs
    header = {}
    for key, val in h5_header.items():
        header[key] = val

    return header

def get_h5_vector(h5_f, key):
    return np.array(h5_f["PartType1/" + key][()])

def set_h5_vector(h5_f, key, val):
    if key in h5_f["PartType1"].keys():
        del h5_f["PartType1/" + key]
    dtype, shape = get_dtype_shape(val)
    h5_f.create_dataset("PartType1/" + key, data=val, dtype=dtype)


def set_h5_header(h5_f, header):
    for key, val in header.items():
        set_h5_header_attr(h5_f, key, val)


def set_h5_header_attr(h5_f, key, val):
    attrs = h5_f["Header"].attrs

    if key in attrs.keys():
        attrs.modify(key, val)
    else:
        dtype, shape = get_dtype_shape(val)
        # if shape == 1 and hasattr(val, "__len__"):
            # val = val[0]
        attrs.create(key, val, shape, dtype=dtype)


def get_dtype_shape(val):
    if hasattr(val, "__len__"):
        return get_arr_dtype_shape(val)
    else:
        shape = 1
        dtype = get_dtype(val)
        return dtype, shape


def get_arr_dtype_shape(arr):
    if isinstance(arr, np.ndarray):
        shape = arr.shape
        if len(shape) == 1:
            shape = shape[0]
        el = arr.flatten()[0]
        dtype = get_dtype(el)
        return dtype, shape

    elif isinstance(arr, (tuple, list)):
        shape = len(arr)
        el = arr[0]
        dtype = get_dtype(el)
        return dtype, shape

    else:
        raise NotImplementedError


def get_dtype(val):
    if isinstance(val, (float, np.floating)):
        return "double"
    if isinstance(val, (int, np.integer)):
        return "int"
    if isinstance(val, (str, np.str_)):
        return "str"
    else:
        raise NotImplementedError(type(val))
