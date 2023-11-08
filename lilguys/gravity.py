import numpy as np
from .units import G

from scipy.optimize import fmin_bfgs
from scipy.special import expi
from scipy.interpolate import interp1d
from dataclasses import dataclass

def numpy_in(f, arg=0):
    def inner(*args, **kwargs):
        args = list(args)
        a = args[arg]
        if isinstance(a, (np.ndarray, tuple, list)):
            args[arg] = np.array(a)
        else:
            raise ValueError("function takes array like arguments")
        return f(*args, **kwargs)

    return inner

@numpy_in
def norm(a, axis=-1):
    return np.sqrt(np.sum(a**2, axis=axis))

@numpy_in
def dist(a, b, axis=-1):
    return norm(a-b, axis=axis)

@numpy_in
def normalize(a, axis=-1):
    return a/norm(a, axis=axis)[:, np.newaxis]

def phi(position, snap):
    """
    Calculates the gravitational potential at r

    Parameters
    ----------
    snap: Snapshot
        A snapshot object describing distribution

    r: vector (3,)
        The position

    epsilon: float
        The gravitational smoothing length

    Returns
    -------
    potential: float
        The gravitational potential 
    """
    r = dist(position, snap.pos)
    r_eff = np.sqrt(r**2 + snap.epsilon**2)
    potential = -G * snap.m * np.sum( 1/r_eff )
    return potential


def grad_phi(position, snap):
    delta_r = position - snap.pos
    r = norm(delta_r)
    r_eff = np.sqrt(r**2 + snap.epsilon**2)
    r_eff = r_eff.reshape(-1, 1)
    return -G * snap.m *np.sum(delta_r / r_eff**3, axis=0)

def min_phi(snap):
    initial_guess = snap.pos[np.argmin(snap.potential)]
    result = fmin_bfgs(phi, initial_guess, fprime=grad_phi, args=(snap,))
    return result


# this section contains monte carlo code

def rand_unit_vector(N=1):
    """
    Creates N random unit vectors
    """
    vect = np.random.normal(0, 1, (N, 3))
    return normalize(vect)

def rho_nfw(r, params):
    """
    The density of an NFW halo
    
    Parameters
    ----------
    r: float
        the radius to evaluate the density at
    rho_0: float
        The characteristic density of the halo
    r_s: float
        The scale radius of the halo
    """
    rho_0 = params.rho_0
    r_s = params.r_s
    return rho_0 / (r/r_s * (1+r/r_s)**2)

def M_nfw(r, params):
    rho_0 = params.rho_0
    r_s = params.r_s
    4*pi*(r_s^2*log(r + r_s) + r_s^3/(r + r_s))*r_s*rho_0




@dataclass
class params_en:
    rho_0: float
    r_s: float
    r_cut: float
    kappa: float = 0.3


def rho_en(r, params):
    rho_0 = params.rho_0
    r_s = params.r_s
    r_cut = params.r_cut
    kappa = params.kappa
    return rho_nfw(r, params) * np.exp(-r/r_cut) / (1+r_s/r_cut)**kappa


def M_en(r, params):
    rho_0 = params.rho_0
    r_s = params.r_s
    r_cut = params.r_cut
    kappa = params.kappa
    pi = np.pi
    return -4*(pi*r*r_cut*r_s**3*np.exp(r/r_cut) + (pi*r_cut*np.exp(r/r_cut) - pi*r_cut)*r_s**4 - (pi*r*r_cut*r_s**3*np.exp(r/r_cut) + pi*r_s**5*np.exp(r/r_cut) + (pi*r + pi*r_cut)*r_s**4*np.exp(r/r_cut))*expi(-(r + r_s)/r_cut)*np.exp(r_s/r_cut) + (pi*r*r_cut*r_s**3*np.exp(r/r_cut) + pi*r_s**5*np.exp(r/r_cut) + (pi*r + pi*r_cut)*r_s**4*np.exp(r/r_cut))*expi(-r_s/r_cut)*np.exp(r_s/r_cut))*rho_0/((r*r_cut*np.exp(r/r_cut) + r_cut*r_s*np.exp(r/r_cut))*((r_cut + r_s)/r_cut)**kappa)

def inverse_M_en(params, N_eval=1000, r_min=None, r_max=None):
    if r_max is None:
        r_max = params.r_cut * 10
    if r_min is None:
        r_min = params.r_s /10
    radii = np.logspace(np.log10(r_min), np.log10(r_max), N_eval)
    mass_profile = M_en(radii, params)
    mass_profile /= mass_profile[-1]
    return interp1d(mass_profile, radii, bounds_error=False, fill_value="extrapolate")

def rand_pos_en(params, Neval, Npoints):
    inverse_M = inverse_M_en(params, Neval)
    rand_M_frac = np.random.uniform(0, 1, Npoints)
    radii = inverse_M(rand_M_frac)[:, np.newaxis]
    return radii * rand_unit_vector(Npoints)


