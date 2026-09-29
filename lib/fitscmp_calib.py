"""Frame calibration: master dark subtraction and correction of isolated hot pixels."""
import numpy as np
from scipy import ndimage

from .fitscmp_io import list_images, load_image

_NEIGHBORS = np.ones((3, 3), bool)
_NEIGHBORS[1, 1] = False  # the 8 neighbors of a pixel


def load_dark(path):
    """Master dark: first 2D image of the file. Returns (array, label)."""
    refs = list_images(path)
    if not refs:
        raise ValueError(f"No 2D image in the dark file {path}")
    return load_image(refs[0]), refs[0].label()


def _noise_model(level, r, nbins=20):
    """Coefficients (a, b) of var(r) = a + b * level (read + shot noise), from the MAD of r in
    quantile bins of the local level; no camera gain is needed."""
    edges = np.unique(np.quantile(level, np.linspace(0.0, 1.0, nbins + 1)))
    lev, var = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        sel = (level >= lo) & (level <= hi)
        if sel.sum() >= 50:
            rs = r[sel]
            lev.append(np.median(level[sel]))
            var.append((1.4826 * np.median(np.abs(rs - np.median(rs)))) ** 2)
    if len(lev) < 2 or np.ptp(lev) == 0:
        return (var[0] if var else np.var(r)), 0.0
    b, a = np.polyfit(lev, var, 1)
    return max(a, 1e-12), max(b, 0.0)


def hot_pixels(img, nsigma=5.0):
    """Boolean mask of isolated hot pixels.

    A pixel is hot if it exceeds the maximum of its 8 neighbors by more than nsigma times the
    local noise, and all neighbors stay below half of its excess over the local background
    (median of the neighbors): the peak of a sampled spot, whose neighbors are comparable,
    is not flagged. The local noise grows with the local level (read + shot noise model
    estimated from the image). Clusters of adjacent hot pixels are not detected.
    """
    finite = np.isfinite(img)
    x = np.where(finite, img, np.median(img[finite]))
    nmax = ndimage.maximum_filter(x, footprint=_NEIGHBORS, mode="nearest")
    nmed = ndimage.median_filter(x, footprint=_NEIGHBORS, mode="nearest")
    r = x - nmed
    level = np.maximum(nmed - np.median(nmed), 0.0)
    a, b = _noise_model(level, r)
    noise = np.sqrt(a + b * level)
    return finite & (x - nmax > nsigma * noise) & (nmax - nmed < 0.5 * r)


def calibrate(img, dark=None, hot_filter=False, nsigma=5.0, saturation=None):
    """Calibrated image, saturation mask and number of corrected hot pixels.

    The saturation mask is computed on the raw values (before dark subtraction), then hot
    pixels are removed from it. Hot pixels are searched after dark subtraction and replaced
    by the median of their 8 neighbors.
    """
    if dark is not None and dark.shape != img.shape:
        raise ValueError(f"Dark shape {dark.shape} differs from image shape {img.shape}")
    sat = None if saturation is None else np.isfinite(img) & (img >= saturation)
    out = img if dark is None else img - dark
    n_hot = 0
    if hot_filter:
        hot = hot_pixels(out, nsigma)
        n_hot = int(hot.sum())
        if n_hot:
            finite = np.isfinite(out)
            med = ndimage.median_filter(np.where(finite, out, np.median(out[finite])),
                                        footprint=_NEIGHBORS, mode="nearest")
            out = np.where(hot, med, out)
            if sat is not None:
                sat &= ~hot  # a saturated hot pixel is not a saturated spot
    return out, sat, n_hot
