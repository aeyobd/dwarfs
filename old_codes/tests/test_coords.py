import pytest
from pytest import approx
import numpy as np

from lilguys import coords

def test_to_galcen():
    gc =  coords.observation(ra = 266.4168166, dec=-29.00782, distance=8.29,
            pm_ra=0, pm_dec=0, radial_velocity=0)
    ph = coords.to_galcen(gc)

    x = ph.x
    y = ph.y
    z = ph.z
    assert x < 1e-2
    assert y < 1e-2
    assert z < 1e-2



def test_to_sky():
    g = coords.phase_point(0,0,0,0,0,0)
    obs = coords.to_sky(g)

    ra = obs.ra
    dec = obs.dec
    assert obs.distance == approx(8.29, rel=1e-4)
    assert obs.ra == approx(266.4168166, rel=1e-4)
    assert obs.dec == approx(-29.00782, rel= 1e-2)


def test_inverse():
    N = 100
    phase = coords.phase_point(
            *np.random.normal(0, 20, (3, N)),
            *np.random.normal(0, 100, (3, N)))

    phase_2 = coords.to_galcen(coords.to_sky(phase))

    assert phase_2.x - phase.x == approx(0, abs=1e-2)
    assert phase_2.y - phase.y == approx(0, abs=1e-2)
    assert phase_2.z - phase.z == approx(0, abs=1e-2)
    assert phase_2.v_z - phase.v_z == approx(0, abs=1e-2)
    assert phase_2.v_y - phase.v_y == approx(0, abs=1e-2)
    assert phase_2.v_x - phase.v_x == approx(0, abs=1e-2)

