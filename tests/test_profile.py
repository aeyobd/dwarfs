import pytest
from pytest import approx

import numpy as np

from lilguys import profile
from lilguys import Snapshot


def test_rho_star_int():
    result = profile.rho_star_int(np.inf, 2, n=4)
    expected = 5_109_350_400 * np.pi
    assert np.isclose(result, expected)

N=10000
pos = np.random.uniform(-1, 1, (N, 3))
vel = np.random.normal(0, 0.05, (N, 3))
pot = np.random.uniform(-1, 1, N)

snap = Snapshot(pos, vel, potential = pot)
x0 = np.random.uniform(3)
v0 = np.random.uniform(3)
shifted = snap.shift(x0, v0)
tol = 2e-2

def test_get_most_bound():
    xb,vb = profile.get_most_bound(shifted)

    assert x0 - xb == approx(0, abs=tol)
    assert v0 - vb == approx(0, abs=tol)

def test_center_snapshot():
    centered = profile.center_snapshot(shifted)
    xb, vb = profile.get_most_bound(centered)
    assert xb == approx(0, abs=2*tol)
    assert vb == approx(0, abs=2*tol)
    

def test_density():
    pass

def test_sort_r():
    snap_sorted = profile.sort_r(snap)
    idx_unsort = np.argsort(np.argsort(snap.r))
    assert np.all(snap.r == snap_sorted.r[idx_unsort])
    assert np.all(snap.IDs == snap_sorted.IDs[idx_unsort])
