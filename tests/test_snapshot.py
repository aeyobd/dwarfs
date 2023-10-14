import pytest
import numpy as np
import os
import h5py

from lilguys import Snapshot
from lilguys import snapshot

Nrows = 100
pos = np.random.uniform(-1, 1, [Nrows, 3])
vel = np.random.uniform(-1, 1, [Nrows, 3])
potentials = np.random.uniform(-1, 0, Nrows)
ids = np.arange(Nrows)
filename = "tests/snapshot_test.hdf5"
newfilename = "tests/snapshot_test_new.hdf5"

def test_get_dtype():
    dtypes = {
            "double": 0.0,
            "int": 123,
            "str": "abc",
            }
    for dtype, val in dtypes.items():
        assert dtype == snapshot.get_dtype(val)

def test_init():
    if os.path.exists(filename):
        os.remove(filename)
    if os.path.exists(newfilename):
        os.remove(newfilename)

def test_save_snapshot():
    snap = Snapshot(pos, vel, ids, potentials)
    assert snap.save(filename)

def test_save_as():
    snap = Snapshot.file(filename)
    assert snap.save(newfilename)
    newsnap = Snapshot.file(newfilename)
    assert np.all(newsnap.pos == pos)

def test_save_all_cols():
    expected_cols = ['Acceleration', 'Coordinates', 'ParticleIDs', 'PertAccel', 'Potential', 'TimeStep', 'Velocities']
    with h5py.File(newfilename) as f:
        cols = list(f["PartType1"].keys())
        assert np.all(np.isin(expected_cols, cols))


def test_save_all_attrs():
    pass


def test_open_snapshot():
    snap = Snapshot.file(filename)
    assert len(snap) == Nrows


def test_xyz():
    snap = Snapshot.file(filename)
    assert np.all(pos[:,0] == snap.x)
    assert np.all(pos[:,1] == snap.y)
    assert np.all(pos[:,2] == snap.z)


def test_transform():
    dp = np.random.normal(0, 1, (3))
    dv = np.random.normal(0, 1, (3))
    dm = np.random.uniform(1, 10)

    snap = Snapshot.file(filename)


