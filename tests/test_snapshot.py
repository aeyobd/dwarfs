import pytest
import numpy as np

from lilguys import Snapshot

Nrows = 100
pos = np.random.uniform(-1, 1, [Nrows, 3])
vel = np.random.uniform(-1, 1, [Nrows, 3])
potentials = np.random.uniform(-1, 0, Nrows)
ids = np.arange(Nrows)
filename = "snapshot_test.hdf5"

def test_save_snapshot():
    snap = Snapshot(pos, vel, ids, potentials)
    assert snap.save("tests/snapshot_test.hdf5")

def test_open_snapshot():
    snap = Snapshot.file("tests/snapshot_000.hdf5")
    assert len(snap.pos)

def test_xyz():
    snap = Snapshot.file("tests/snapshot_000.hdf5")
    pos = snap.pos
    assert np.all(pos[:,0] == snap.x)
    assert np.all(pos[:,1] == snap.y)
    assert np.all(pos[:,2] == snap.z)
