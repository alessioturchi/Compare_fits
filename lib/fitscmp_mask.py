"""Masked-pixel handling for the Fourier analysis: apodized FFT weights and inpainting."""
import warnings

import numpy as np
from astropy.convolution import Gaussian2DKernel, interpolate_replace_nans
from astropy.utils.exceptions import AstropyUserWarning
from scipy.ndimage import distance_transform_edt

NAN_MODES = ("mask", "interp")


def fft_weights(mask, window=True, taper=0.0):
    """Weight map for the FFT: optional 2D Hann window times the apodized validity mask.

    Invalid pixels get weight 0; valid pixels closer than `taper` pixels to an
    invalid one are down-weighted with a cosine ramp going from 0 to 1.
    """
    ny, nx = mask.shape
    w = np.outer(np.hanning(ny), np.hanning(nx)) if window else np.ones(mask.shape)
    if taper > 0 and not mask.all():
        dist = distance_transform_edt(mask)  # distance to the nearest invalid pixel (0 on invalid)
        w = w * np.where(dist < taper, 0.5 * (1.0 - np.cos(np.pi * dist / taper)), 1.0)
    return w * mask


def inpaint(img, sigma=1.0):
    """Fill NaN pixels by Gaussian-kernel interpolation; holes larger than the kernel stay NaN."""
    if np.isfinite(img).all():
        return img
    with warnings.catch_warnings():  # residual NaNs in large holes are expected here
        warnings.simplefilter("ignore", AstropyUserWarning)
        return interpolate_replace_nans(img, Gaussian2DKernel(x_stddev=sigma), boundary="extend")


def inpaint_common(a, b, sigma=1.0):
    """Fill pixels invalid in A or B by interpolation, identically in both images."""
    bad = ~(np.isfinite(a) & np.isfinite(b))
    if not bad.any():
        return a, b
    return tuple(inpaint(np.where(bad, np.nan, x), sigma) for x in (a, b))
