from . import Snapshot
import numpy as np
import pandas as pd
from glob import glob
from os import path

class Output:
    def __init__(self, directory):
        filenames = glob(path.join(directory + "/*.hdf5"))
        self.Nt = len(filenames)
        snap0 = Snapshot.file(filenames[0])
        self.Np = len(snap0)

        self.pos = np.empty((self.Np, self.Nt, 3))
        self.vel = np.empty((self.Np, self.Nt, 3))
        self.t = np.empty(self.Nt)
        self.pot = np.empty((self.Np, self.Nt))
        self.m = snap0.m

        for i in range(self.Nt):
            snap = Snapshot.file(filenames[i])
            for j in range(self.Np):
                idx =  snap.IDs[j]
                self.pos[idx, i, :] = snap.pos[j]
                self.vel[idx, i, :] = snap.vel[j]
                self.pot[idx, i] = snap.potential[j]
            self.t[i] = snap.header["Time"]
