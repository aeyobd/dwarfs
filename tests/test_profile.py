import pytest
from pytest import approx

import numpy as np
from scipy.stats import chi2

from lilguys import profile
from lilguys import Snapshot
from lilguys import units



def test_rho_star_int():
    result = profile.rho_star_int(np.inf, 2, n=4)
    expected = 5_109_350_400 * np.pi
    assert np.isclose(result, expected)


class TestSnapshotUtils:
    N=10000
    pos = np.random.uniform(-1, 1, (N, 3))
    vel = np.random.normal(0, 0.05, (N, 3))
    m=1

    snap = Snapshot(pos, vel, m=m, epsilon=1/np.sqrt(N))
    snap.calc_potential()
    x0 = np.random.normal(0, 1, (3))
    v0 = np.random.uniform(0, 1, (3))
    shifted = snap.shift(x0, v0)

    dx = np.sqrt(np.mean(snap.r**2)) / np.sqrt(N)

    def test_get_center(self):
        xb, vb = profile.get_center(self.shifted)
        assert self.x0 - xb == approx(0, abs=5*self.dx)
        assert self.v0 - vb == approx(0, abs=5*self.dx)

    def test_center_snapshot(self):
        centered = profile.center_snapshot(self.shifted)
        xb, vb = profile.get_center(centered)
        assert xb == approx(0, abs=5*self.dx)
        assert vb == approx(0, abs=5*self.dx)
    

    def test_sort_r(self):
        sort_snap = profile.sort_r(self.snap)
        idx_unsort = np.argsort(np.argsort(self.snap.r))
        assert np.all(sort_snap.r[1:] > sort_snap.r[:-1])
        assert np.all(self.snap.IDs == sort_snap.IDs[idx_unsort])

    def test_get_KE(self):
        KE = profile.get_KE(self.snap)
        E_exp = 1/2 * self.snap.v**2
        assert  KE - E_exp == approx(0, abs=1e-8)




def test_V_circ():
    N=10000
    Rmax = 10
    x = np.random.uniform(0, Rmax, N)
    y = np.random.uniform(0, Rmax, N)
    r = np.sort(np.sqrt(x**2 + y**2))
    r = r[r<Rmax]
    N = len(r)
    M = np.arange(N) + 0.5
    Vcirc_exp = np.sqrt(units.G * N * r)/Rmax
    Vcirc = profile.V_circ(M, r)
    Vcirc_err = np.abs(Vcirc)/np.sqrt(1+np.arange(N))
    # these ar eorder 1
    resid = (Vcirc - Vcirc_exp)
    chi_sq = np.sum(np.square(resid/Vcirc_err))
    assert chi_sq < chi2.ppf(0.99, N)


def test_get_L():
    x = np.array([(1,0,0), (1,0,0), (3,0,0)])
    v = np.array([(0,0,0), (0,2,0), (0,0,1)])
    snap = Snapshot(x, v, m=1)
    L = profile.get_L(snap)
    L_exp = np.array([
        (0,0,0),
        (0,0,2),
        (0,-3,0)
        ])
    assert L - L_exp  == approx(0, abs=1e-8)

def normalize(a, axis=-1):
    l2 = np.linalg.norm(a, axis=axis)
    l2[l2==0] = 1
    return a / np.expand_dims(l2, axis)

def rand_unit(shape=(1, 3)):
    rand_vec = np.random.normal(0, 1, shape)
    return normalize(rand_vec)

def z_score(actual, expected, error):
    return np.abs(actual - expected) / error

def norm(vec):
    return np.sum(np.square(vec), axis=-1)

class TestUniformProfile:
    N = 10_000
    Rmax = 1
    pos = np.random.uniform(-Rmax, Rmax, (N, 3))
    vel = pos
    m=1
    Nbins = 10
    snap = Snapshot(pos, vel, m=m, epsilon=1/np.sqrt(N))

    snap.calc_potential()
    profile.center_snapshot(snap, inplace=True)
    
    filt = snap.r < Rmax
    snap.filter(filt, inplace=True)
    N = len(snap)

    prof = profile.Profile(snap, r_bins=Nbins, E_bins = Nbins, center=False)
    V = 4/3*np.pi * Rmax**3
    rho = m * N / V

    def test_create_prof(self):
        assert self.prof is not None
        assert len(self.prof.r) == self.Nbins
        assert len(self.prof.nu_DM) == self.Nbins
        assert len(self.prof.M) == self.Nbins
        assert self.prof.m == self.m

    def test_r(self):
        assert np.all(self.prof.r[1:] > self.prof.r[:-1])
        assert self.prof.r[-1] <= np.sqrt(3) * self.Rmax

    def test_rho(self):
        rho_exp = self.rho
        rho = self.prof.nu_DM
        rho_err = self.prof.nu_DM_err
        print(rho)
        print(rho_exp)
        print(np.max(self.snap.r))
        print(self.prof.r_bins)
        z = z_score(rho_exp, rho, rho_err)

        assert np.nanmean(z**2) < 2
        assert np.all(np.nan_to_num(rho_err, 0) <= rho)

    def test_M(self):
        M_exp = self.m * self.N * (self.prof.r/self.Rmax)**3
        M = self.prof.M
        M_err = self.prof.M_err
        
        print(M)
        print(M_exp)
        z = z_score(M_exp, M, M_err)
        assert np.all(z < 5)
        assert np.all(M_err <= M)

    def test_V_circ(self):
        V_0 = np.sqrt(self.m * self.N / self.Rmax)
        V_circ = self.prof.V_circ
        V_circ_err = self.prof.V_circ_err
        V_circ_exp = self.prof.r/self.Rmax * V_0

        print(V_circ)
        print(V_circ_exp)
        z = z_score(V_circ_exp, V_circ, V_circ_err)
        assert np.all(z < 5)



