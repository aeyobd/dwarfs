import numpy as np
from scipy.spatial import KDTree
import matplotlib.pyplot as plt
import matplotlib as mpl


def mag_image(x, y, im_kwargs={}, **kwargs):
    image, x_edges, y_edges = kernel_hist2d(x, y, **kwargs)

    imin = np.min(image[image>0])
    imin = max(imin, 0.05)
    norm = mpl.colors.LogNorm(imin, np.max(image))
    extent=[min(x_edges), max(x_edges), min(y_edges), max(y_edges)]

    plt.imshow(image, origin="lower", cmap="gray_r",
            extent=extent, norm=norm, **im_kwargs)


def kernel_hist2d(x, y, knn=10, depth=5, eta=2, **kwargs):
    """
    Creates a kernel-smoothed histogram from the aray of points

    Parameters
    ----------
    x:  `np.ndarray`
        The array of x values
    y:  `np.ndaray`
        The array of y values
    knn:    `int`
        The number of neighbors to use for density estimation
    eta:    `float`
        The kernel width scaling factor
    depth:  `float`
        Maximum number of points to smear

    Passed to `np.histogram2d`:
    bins : `int` [default = 10]
    range: array like, 2x2
    """

    X = np.array([x, y]).transpose()
    tree = KDTree(X)

    image, x_edges, y_edges = np.histogram2d(X[:,0], X[:,1], **kwargs)
    x_mid = bin_mids(x_edges)
    y_mid = bin_mids(y_edges)

    image_old = image.copy()


    dx = x_edges[1] - x_edges[0]

    coords = np.where((image > 0))
    itercoords = np.transpose(coords)
    weights = np.minimum(image[coords], depth)
    image[coords] = image[coords] - weights

    for i, j, w in zip(*coords, weights):
        dist = tree.query((x_mid[i], y_mid[j]), k=knn)[0][-1]
        sigma = eta * dist/dx / np.sqrt(knn)
        kern = gaussian_kernel(sigma)
        add_overlay(image, w*kern, i, j)

    return image, x_edges, y_edges


def gaussian_kernel(sigma=1, side_length=None):
    if side_length is None:
        half_length = int(np.ceil(5*sigma))
    else:
        half_length = (side_length - 1)/2

    if sigma==0:
        return np.zeros(side_length)

    r = np.arange(-half_length, half_length+1, 1)
    gauss = np.exp(-1/2 * np.square(r)/np.square(sigma))
    kernel = np.outer(gauss, gauss)
    return kernel / np.sum(kernel)


def add_overlay(image, overlay, center_x, center_y, inplace=False):
    """
    Add a subregion from a smaller array to a given position in a full image.

    Parameters:
        image (numpy.ndarray): The full image where the subregion will be added.
        overlay (numpy.ndarray): The smaller array to be added.
        center_x (int): X-coordinate of the center of the subregion in the full image.
        center_y (int): Y-coordinate of the center of the subregion in the full image.

    Returns:
        numpy.ndarray: The full image with the subregion added.

    This function adds a subregion from the smaller_array to the full_image. It centers the subregion on the specified
    pixel coordinates (center_x, center_y) within the full image while ignoring parts that extend beyond the
    boundaries of the full image.

    Example usage:
    full_image = np.zeros((20, 20))
    smaller_array = np.ones((5, 5))
    result_image = add_subregion_to_image(full_image, smaller_array, 10, 10)
    """
    if not inplace:
        image = np.copy(image)
    side_length = overlay.shape[0]
    half_length = side_length // 2

    start_x = max(0, center_x - half_length)
    end_x = min(image.shape[0], center_x + half_length + 1)
    start_y = max(0, center_y - half_length)
    end_y = min(image.shape[1], center_y + half_length + 1)

    subregion = image[start_x:end_x, start_y:end_y]

    overlay_start_x = max(0, half_length - center_x)
    overlay_end_x = overlay_start_x + (end_x - start_x)
    overlay_start_y = max(0, half_length - center_y)
    overlay_end_y = overlay_start_y + (end_y - start_y)

    subregion += overlay[overlay_start_x:overlay_end_x, overlay_start_y:overlay_end_y]

    image[start_x:end_x, start_y:end_y] = subregion

    return image

def bin_mids(bins):
    return (bins[1:] + bins[:-1])/2
