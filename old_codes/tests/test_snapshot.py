import pytest
import numpy as np
import os
import h5py

from lilguys import Snapshot
from lilguys import snapshot


filename = "tests/snapshot_test.hdf5"
newfilename = "tests/snapshot_test_new.hdf5"

@pytest.fixture
def clean_dir():
    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(newfilename):
        os.remove(newfilename)

class TestSnapshot():
    Nrows = 100
    pos = np.random.uniform(-1, 1, [Nrows, 3])
    vel = np.random.uniform(-1, 1, [Nrows, 3])
    potentials = np.random.uniform(-1, 0, Nrows)
    ids = np.arange(Nrows)
    filename = filename
    newfilename = newfilename

    def test_save_snapshot(self):
        snap = Snapshot(self.pos, self.vel, self.ids, self.potentials)
        assert snap.save(self.filename)

    @pytest.fixture
    def snap(self):
        return Snapshot.file(self.filename)

    def test_save_as(self, snap):
        assert snap.save(self.newfilename)
        newsnap = Snapshot.file(self.newfilename)
        assert np.all(newsnap.pos == self.pos)

    def test_copy(self, snap):
        snap1 = snap.copy()
        dx = np.random.uniform(3)
        snap1.shift(dx, inplace=True)
        assert np.all(snap1.pos == snap.pos+dx)

    def test_shift(self, snap):
        dx = np.random.uniform(3)
        dv = np.random.uniform(3)
        shifted = snap.shift(dx, dv)
        assert np.all(snap.pos + dx == shifted.pos)
        assert np.all(snap.vel + dv == shifted.vel)


    def test_save_all_cols(self, snap):
        expected_cols = ['Acceleration', 'Coordinates', 'ParticleIDs', 'Potential', 'Velocities']
        with h5py.File(newfilename) as f:
            cols = list(f["PartType1"].keys())
            assert np.all(np.isin(expected_cols, cols))


    def test_save_all_attrs(self):
        pass


    def test_open_snapshot(self, snap):
        assert len(snap) == self.Nrows


    def test_xyz(self, snap):
        assert np.all(self.pos[:,0] == snap.x)
        assert np.all(self.pos[:,1] == snap.y)
        assert np.all(self.pos[:,2] == snap.z)

    def test_vxyz(self, snap):
        assert np.all(self.vel[:,0] == snap.v_x)
        assert np.all(self.vel[:,1] == snap.v_y)
        assert np.all(self.vel[:,2] == snap.v_z)
