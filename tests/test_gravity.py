import pytest
from pytest import approx

import numpy as np

from lilguys import gravity
from lilguys.units import G
import lilguys


def test_norm():
    assert gravity.norm([0, 3, 4]) == 5.
    assert gravity.norm(np.zeros(3)) == 0
    assert gravity.norm(np.ones(3)) == np.sqrt(3)
    assert np.all(gravity.norm(np.random.normal(0, 1, (100, 3))) >= 0)

def test_dist():
    N = 100
    a = np.random.normal(0, 1, (N, 3))
    b = np.random.normal(0, 1, (N, 3))
    c = np.random.normal(0, 1, (N, 3))

    assert gravity.dist((0,0,0), (0,3,4)) == 5
    assert np.all(gravity.dist(a, a) == 0)
    assert np.all(gravity.dist(a, b) == gravity.dist(b, a))
    assert np.all(gravity.dist(a, b) >= 0)
    assert np.all(gravity.dist(a, b) + gravity.dist(b, c) >= gravity.dist(a, c))

@pytest.fixture
def twobody():
    pos = np.array([[0,0,0], [1,0,0]])
    vel = np.random.normal(0, 1, (2,3))
    m = 1
    epsilon = 0.01

    snap = lilguys.Snapshot(pos, vel, m=m)
    snap.epsilon = epsilon
    return snap

def test_phi(twobody):
    p = [0.5, 0, 0]
    phi_exp = - 2 * G * twobody.m / 0.5

    assert gravity.phi(p, twobody) == approx(phi_exp, rel=1e-3)
    assert gravity.grad_phi(p, twobody) == approx(0, abs=1e-3)


    p = [0, 0, 0]
    phi_exp = - G * twobody.m * (1 + 1/twobody.epsilon)
    assert gravity.phi(p, twobody) == approx(phi_exp, rel=1e-3)
    assert gravity.grad_phi(p, twobody) == approx([G*twobody.m/1**2, 0, 0], abs=1e-2)

@pytest.fixture
def nbody():
    N = 300
    M = 30
    pos = np.random.normal(1, 0.1, (N, 3))
    pos = np.concatenate([pos, np.random.uniform(-20, 20, (M, 3))])
    vel = np.random.normal(0, 1, (N + M, 3))
    m = 1
    epsilon = 0.1/np.sqrt(N)

    snap = lilguys.Snapshot(pos, vel, m=m)
    snap.epsilon = epsilon
    return snap

def test_minphi(nbody):
    p0 = np.ones(3)
    p = gravity.min_phi(nbody)

    sigma = 10*0.1/np.sqrt(len(nbody)) 
    assert p == approx(p0, abs=5*sigma)

