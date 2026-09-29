"""Image registration: sub-pixel translation estimate and Fourier-shift resampling."""
import math

import numpy as np
from scipy.ndimage import shift as ndi_shift
from skimage.registration import phase_cross_correlation

from .fitscmp_mask import fft_weights, inpaint, inpaint_common


def estimate_shift(a, b, upsample=100, taper=16.0):
    """Shift (dy, dx) [pix] to apply to B to register it onto A.

    Cross-correlation peak refined by upsampled DFT (Guizar-Sicairos et al. 2008), on
    mean-subtracted images multiplied by the common apodized mask and a Hann window.
    Plain (not phase-normalized) correlation is used: phase normalization weights
    noise-dominated frequencies equally and degrades band-limited images.
    """
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")
    a, b = inpaint_common(a, b)
    mask = np.isfinite(a) & np.isfinite(b)
    w = fft_weights(mask, True, taper)
    if not np.any(w > 0):
        raise ValueError("No valid pixel for registration")
    xa, xb = (np.where(mask, x - np.mean(x[mask]), 0.0) * w for x in (a, b))
    shift, _, _ = phase_cross_correlation(xa, xb, upsample_factor=int(upsample), normalization=None)
    return float(shift[0]), float(shift[1])


def apply_shift(img, dy, dx, pad=16):
    """Shift img by (dy, dx) pixels (content moves towards +y, +x) with the Fourier shift theorem.

    The image is mirror-padded so that the circular shift does not wrap the opposite
    border into the field. Pixels entering from outside the field, and pixels whose
    value depends on NaN input pixels, are returned as NaN.
    """
    ny, nx = img.shape
    if abs(dy) >= ny or abs(dx) >= nx:
        raise ValueError(f"Shift ({dy}, {dx}) exceeds the image size {img.shape}")
    if dy == 0 and dx == 0:
        return img
    valid = np.isfinite(img)
    x = img
    if not valid.all():  # the FFT needs finite values: inpaint small holes, mean elsewhere
        x = inpaint(img)
        x = np.where(np.isfinite(x), x, np.mean(img[valid]))
    py = min(math.ceil(abs(dy)) + pad, ny - 1)
    px = min(math.ceil(abs(dx)) + pad, nx - 1)
    x = np.pad(x, ((py, py), (px, px)), mode="reflect")
    fy = np.fft.fftfreq(x.shape[0])[:, None]
    fx = np.fft.fftfreq(x.shape[1])[None, :]
    x = np.real(np.fft.ifft2(np.fft.fft2(x) * np.exp(-2j * np.pi * (fy * dy + fx * dx))))
    x = x[py:py + ny, px:px + nx]

    if not valid.all():  # shifted NaN mask: any contribution from an invalid pixel
        x[ndi_shift((~valid).astype(float), (dy, dx), order=1, mode="constant") > 0] = np.nan
    iy, ix = math.ceil(abs(dy)), math.ceil(abs(dx))  # strips with no source pixels
    if dy > 0:
        x[:iy] = np.nan
    elif dy < 0:
        x[ny - iy:] = np.nan
    if dx > 0:
        x[:, :ix] = np.nan
    elif dx < 0:
        x[:, nx - ix:] = np.nan
    return x
