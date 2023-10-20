import h5py
import numpy as np
from os import path

from .hdfutils import get_h5_vector, get_h5_header, set_h5_vector, set_h5_header, make_default_header


class Snapshot:
    _filename = None
    gadget4 = False

    def __init__(self, positions, velocities, 
            IDs=None, potential=None, accelerations=None,
            header = {}, m=None):
        self.pos = positions
        self.vel = velocities
        if IDs is None:
            IDs = np.arange(len(self))
        self.acc = accelerations
        self.IDs = IDs
        if potential is None:
            potential = np.zeros(len(self))
        self.potential = potential
        if header == {}:
            N = len(self)
            header = make_default_header(N, m)
        self.header = header

    @classmethod
    def file(cls, filename):
        with h5py.File(filename, "r") as f:
            keys = f["PartType1"].keys()
            header = get_h5_header(f)

            pos = get_h5_vector(f, "Coordinates")
            vel = get_h5_vector(f, "Velocities")
            IDs = get_h5_vector(f, "ParticleIDs")

            if "Acceleration" in keys:
                acc = get_h5_vector(f, "Acceleration")
            else:
                acc = None
            if "Potential" in keys:
                pot = get_h5_vector(f, "Potential")
            else:
                pot = None

        c = cls(pos, vel, IDs=IDs, header=header, potential=pot, 
                accelerations=acc)
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
                set_h5_vector(f, "Acceleration", np.zeros(len(self)))
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




