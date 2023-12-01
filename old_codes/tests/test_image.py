import pytest
from pytest import approx

from lilguys import image
import numpy as np


def test_gaussian_kern_norm():
    s = np.random.uniform(1, 5, 100)
    area = [np.sum(image.gaussian_kernel(x)) for x in s]

    assert np.array(area) == approx(1, abs=1e-5)

def test_gaussian_kernel_symmetry():
    sigma = 2
    s = 5*sigma
    kern = image.gaussian_kernel(sigma)
    l = 2*s + 1
    assert kern.shape == (l, l)

    for i in range(0, s):
        for j in range(0, s):
            assert kern[i, j] == kern[j, i]
            assert kern[i, j] == kern[l-1 - i, l-1 -j]

def test_gaussian_kernel_profile():
    from scipy.stats import norm
    sigma = 1
    side_length = 21
    kern = image.gaussian_kernel(sigma, side_length=side_length)
    middle = (side_length-1)//2 
    prof = kern[middle]
    prof /= np.sum(prof)

    xs = np.arange(-middle, middle+1)

    std = np.sqrt(np.sum(xs**2 * prof))
    mean = np.sum(xs * prof)
    lim = prof[0] + prof[-1]
    assert std == approx(1, abs=1e-2)
    assert lim == approx(0, abs=1e-5)
    assert mean == approx(0, abs=1e-5)


def test_overlay_image():
    N = 10
    img = np.ones((N,N))
    overlay = np.arange(9).reshape(3,3)
    result = image.add_overlay(img, overlay, 0, 0)
    assert result.shape == (N, N)
    assert np.sum(result - img) == 4+5+7+8


def test_overlay_image_boundaries():
    pass

def test_kern_hist2d_norm():
    pass
def test_kern_hist2d_accuracy():
    pass
