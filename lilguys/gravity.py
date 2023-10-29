import numpy as np
from .units import G

from scipy.optimize import fmin_bfgs


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
    return a/norm(a, axis=axis)

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


