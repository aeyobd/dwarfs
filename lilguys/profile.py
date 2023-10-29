import numpy as np
from scipy.integrate import quad
from scipy.spatial import KDTree
from .units import G
from .gravity import min_phi


def rho_star(r, r_scale, n):
    return np.exp(-(r/r_scale)**(1/n))

def rho_star_int(r_max, r_scale, n):
    return 4*np.pi * quad(lambda r: r**2 * rho_star(r, r_scale, n), 0, r_max)[0]


def V_circ(M, r):
    return np.where(r>0, np.sqrt(G*M/r), 0)

def most_bound(snap, percentile=0.2):
    E = snap.potential
    filt = E < np.percentile(E, percentile)
    return snap.filter(filt)

def get_center(snap, eta=5):
    p0 = snap.pos[np.argmin(snap.potential)]
    s1 = snap.shift(-p0)
    sigma = eta*snap.epsilon
    weights = np.exp(-s1.r**2 / sigma**2)
    v0 = np.average(snap.vel, weights=weights, axis=0)
    return p0, v0


def center_snapshot(snap, inplace=False, verbose=False):
    p0, v0 = get_center(snap)
    if verbose:
        print(f"shifting by {p0}, {v0}")
    return snap.shift(-p0, -v0, inplace=inplace)


def get_KE(snap):
    E_kin = 0.5*snap.v**2
    return E_kin

def get_E_local(snap):
    return get_KE(snap) + snap.potential

def get_Etot(snap):
    return np.sum(0.5*snap.potential + snap.ext_potential + get_KE(snap), axis=-1)

def get_Vmax(snap):
    r = sort_r(snap).r
    M = np.arange(len(snap)) * snap.m
    return np.max(V_circ(M, r))


def get_L(snap, p0=np.zeros(3), v0 = np.zeros(3)):
    return np.cross((snap.pos-p0), (snap.vel - v0))

def sort_r(snap):
    idx = np.argsort(snap.r)
    snap_sorted = snap.copy()
    snap_sorted.pos = snap.pos[idx]
    snap_sorted.vel = snap.vel[idx]
    snap_sorted.IDs = snap.IDs[idx]
    snap_sorted.potential = snap.potential[idx]
    snap_sorted.ext_potential = snap.ext_potential[idx]

    return snap_sorted


class Profile:
    def __init__(self, snap, r_bins=20, E_bins=20, eps_r = 1e-3, center=True):
        if center:
            self.snap = sort_r(center_snapshot(snap))
        else:
            self.snap = sort_r(snap)
        filt = self.snap.r > eps_r
        self.snap = self.snap.filter(filt)

        self.m = snap.m
        self.create_r_bins(r_bins)

        self.compute_masses()
        self.compute_V_circ()
        self.compute_density()
        self.psi = np.interp(self.r, self.snap.r, self.snap.potential)

    def create_r_bins(self, Nbins):
        log_r_min = np.log10(np.min(self.snap.r))
        log_r_max = np.log10(np.max(self.snap.r))
        self.r_bins = np.logspace(log_r_min, log_r_max, Nbins+1)
        self.r = 0.5 * (self.r_bins[1:] + self.r_bins[:-1])


    def compute_V_circ(self):
        self.V_circ = np.sqrt(G * self.M / self.r)
        self.V_circ_err = self.M_err / self.M * self.V_circ
        return self.V_circ

    def compute_masses(self):
        N = len(self.snap.r)
        snap_M_r = self.m * np.arange(N)
        snap_M_err = self.m * np.sqrt(np.arange(N))

        self.M = np.interp(self.r, self.snap.r, snap_M_r)
        self.M_err = np.interp(self.r, self.snap.r, snap_M_err)
        

    def compute_density(self):
        DM_counts, _ = np.histogram(self.snap.r, bins=self.r_bins)
        self.dV = 4/3 * np.pi * (self.r_bins[1:]**3 - self.r_bins[:-1]**3)
        self.nu_DM = DM_counts * self.m / self.dV
        self.nu_DM_err = np.sqrt(DM_counts) * self.m / self.dV
        self.nu_DM_err[DM_counts == 0] = np.nan


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
        self.E = np.linspace(E_min, E_max, Nbins + 1)

    def __len__(self):
        return len(self.r)

