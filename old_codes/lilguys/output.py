from . import Snapshot
import numpy as np
import pandas as pd
from glob import glob
from os import path
from . import hdfutils

class Output:
    def __init__(self, directory):
        filenames = glob(path.join(directory + "/*.hdf5"))
        self.filenames = filenames
        self.Nt = len(filenames)
        snap0 = Snapshot.file(filenames[0])
        self.Np = len(snap0)

        self.pos = np.empty((self.Nt, self.Np, 3))
        self.vel = np.empty((self.Nt, self.Np, 3))
        self.t = np.empty(self.Nt)
        self.potential = np.empty((self.Nt, self.Np))
        self.ext_potential = np.empty((self.Nt, self.Np))
        self.m = snap0.m
        self.IDs = snap0.IDs
        self.epsilon = hdfutils.get_epsilon(directory)

        for i in range(self.Nt):
            snap = Snapshot.file(filenames[i])
            for j in range(self.Np):
                ID =snap.IDs[j]
                idx = np.argwhere(self.IDs == ID)[0][0]
                self.pos[i, idx, :] = snap.pos[j]
                self.vel[i, idx, :] = snap.vel[j]
                self.potential[i, idx] = snap.potential[j]
                self.ext_potential[i, idx] = snap.ext_potential[j]
            self.t[i] = snap.header["Time"]

    @property
    def v(self):
        return np.sqrt(np.sum(self.vel**2, axis=-1))

    @property
    def r(self):
        return np.sqrt(np.sum(self.pos**2, axis=-1))

    def __getitem__(self, idx):
        snap =  Snapshot.file(self.filenames[idx])
        snap.epsilon = self.epsilon
        return snap
