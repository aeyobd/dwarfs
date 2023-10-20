import numpy as np
from scipy.integrate import quad
from scipy.spatial import KDTree
from .units import G


def rho_star(r, r_scale, n):
    return np.exp(-(r/r_scale)**(1/n))

def rho_star_int(r_max, r_scale, n):
    return 4*np.pi * quad(lambda r: r**2 * rho_star(r, r_scale, n), 0, r_max)[0]

def get_most_bound(snap, knn=10):
    X_tree = KDTree(snap.pos)
    k_d, k_i = X_tree.query(snap.pos, k=knn)
    idx = np.argmin(np.mean(k_d, axis=1))
    return snap.pos[idx], snap.vel[idx]

def get_most_bound_old(snap, min_bound=0, verbose=False):
    if snap.potential is None or np.min(snap.potential) >= 0:
        p0 = np.mean(snap.pos, axis=0)
        v0 = np.mean(snap.vel, axis=0)
        return p0, v0

    idx = np.argmin(snap.potential)
    p0 = snap.pos[idx, :]
    v0 = snap.vel[idx, :]


    filt = snap.potential + 0.5 * np.sum((snap.vel-v0)**2, axis=1) < 0
    if verbose:
        print("bound fraction: ", np.mean(filt))

    if np.sum(filt) > min_bound:
        p0 = np.mean(snap.pos[filt], axis=0)
        v0 = np.mean(snap.vel[filt], axis=0)
    return p0, v0

def center_snapshot(snap, inplace=False, verbose=True):
    p0, v0 = get_most_bound(snap)
    if verbose:
        print(f"shifting by {p0}, {v0}")
    return snap.shift(-p0, -v0, inplace=inplace)

def get_energy(snap):
    E_kin = snap.m * snap.v**2
    Etot = E_kin + snap.potential
    return Etot

def sort_r(snap):
    idx = np.argsort(snap.r)
    snap_sorted = snap.copy()
    snap_sorted.pos = snap.pos[idx]
    snap_sorted.vel = snap.vel[idx]
    snap_sorted.IDs = snap.IDs[idx]
    snap_sorted.potential = snap.potential[idx]


    return snap_sorted

class Profile:
    def __init__(self, snap, r_bins=20, E_bins=20, eps_r = 1e-3):
        self.snap = sort_r(center_snapshot(snap))
        self.snap_filt = self.snap.r > eps_r

        self.m = snap.m
        self.create_r_bins(r_bins)

        self.compute_masses()

        self.compute_density()
        self.psi = np.interp(self.r, self.snap.r[self.snap_filt], 
                self.snap.potential[self.snap_filt])

    def create_r_bins(self, Nbins):
        log_r_min = np.log10(np.min(self.snap.r[self.snap_filt]))
        log_r_max = np.log10(np.max(self.snap.r[self.snap_filt]))
        self.r_bins = np.logspace(log_r_min, log_r_max, num=Nbins)
        self.r = 0.5 * (self.r_bins[1:] + self.r_bins[:-1])



    def compute_masses(self):
        N = len(self.snap.r[self.snap_filt])
        snap_M_r = self.m * np.arange(N)
        self.M = np.interp(self.r, self.snap.r[self.snap_filt], snap_M_r)
        

    def compute_density(self):
        DM_counts, _ = np.histogram(self.snap.r[self.snap_filt], bins=self.r_bins)
        self.dV = 4/3 * np.pi * (self.r_bins[1:]**3 - self.r_bins[:-1]**3)
        self.nu_DM = DM_counts * self.m / self.dV

    def set_stars(self, IDs):
        pass
        

    def add_stars(self, n, r_scale):
        self.M_star_tot = rho_star_int(np.inf, r_scale_star, n=n)
        self.M_star = np.array([
            rho_star_int(r, r_scale_star, n) for r in self.r])
        self.nu_star_assumed = rho_star(self.r, r_scale, n)
        self.r_scale = r_scale
        self.n = n

    def create_E_bins(self, Nbins):
        E_max = self.psi[0]
        E_min = E_max/Nbins
        self.E = np.linspace(E_min, E_max, Nbins)

