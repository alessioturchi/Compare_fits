"""Numerical comparison of two 2D images: direct subtraction and Fourier analysis."""
import numpy as np

from .fitscmp_mask import NAN_MODES, fft_weights, inpaint_common

NORMALIZATIONS = ("none", "mean", "median", "max", "sum")
_NORM_FUNC = {"mean": np.nanmean, "median": np.nanmedian, "max": np.nanmax, "sum": np.nansum}


def normalize(img, mode="none", mask=None):
    """Divide the image by a scalar statistic (flux normalization).

    If mask is given, the statistic is computed only on img[mask].
    """
    if mode == "none":
        return img
    if mode not in _NORM_FUNC:
        raise ValueError(f"Unknown normalization '{mode}'")
    f = _NORM_FUNC[mode](img if mask is None else img[mask])
    if not np.isfinite(f) or f == 0:
        raise ValueError(f"Normalization '{mode}' gives invalid factor {f}")
    return img / f


def image_stats(img):
    """Basic statistics of the finite pixels of an image, in its original units."""
    v = img[np.isfinite(img)]
    if v.size == 0:
        return {"n_valid": 0, "n_nan": int(img.size)}
    return {"min": v.min(), "max": v.max(), "mean": v.mean(), "median": np.median(v),
            "std": v.std(), "sum": v.sum(), "n_valid": int(v.size), "n_nan": int(img.size - v.size)}


def difference(a, b, mode="none"):
    """Return normalized images, their difference A-B and statistics on the common finite mask."""
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")
    ok = np.isfinite(a) & np.isfinite(b)  # common mask
    if not ok.any():
        raise ValueError("No pixel is finite in both images")
    a, b = normalize(a, mode, ok), normalize(b, mode, ok)  # same pixel set for both factors
    d = a - b
    rms_d = np.sqrt(np.mean(d[ok] ** 2))
    with np.errstate(divide="ignore", invalid="ignore"):  # constant images give NaN/inf
        stats = {
            "mean": np.mean(d[ok]),
            "std": np.std(d[ok]),
            "rms": rms_d,
            "max_abs": np.max(np.abs(d[ok])),
            "rms_rel": rms_d / np.std(a[ok]),  # RMS(A-B) / std(A)
            "pearson": np.corrcoef(a[ok], b[ok])[0, 1],
            "n_valid": int(ok.sum()),
            "masked": 1.0 - ok.mean(),  # fraction of pixels not finite in both images
        }
    return a, b, d, stats


def _log10(x):
    """log10 with NaN for non-positive values."""
    with np.errstate(divide="ignore", invalid="ignore"):
        y = np.log10(x)
    y[~np.isfinite(y)] = np.nan
    return y


def _centered_fft(img, w):
    """Centered 2D FFT of the mean-subtracted image multiplied by the weight map w.

    The mean is computed on pixels with w > 0; pixels with w == 0 are ignored.
    """
    valid = w > 0
    x = np.where(valid, img - np.mean(img[valid]), 0.0)
    return np.fft.fftshift(np.fft.fft2(x * w))


def _freq_grid(shape):
    """Centered frequency axes (cycles/pixel) and radial frequency map."""
    fy = np.fft.fftshift(np.fft.fftfreq(shape[0]))
    fx = np.fft.fftshift(np.fft.fftfreq(shape[1]))
    fxx, fyy = np.meshgrid(fx, fy)
    return fx, fy, np.hypot(fxx, fyy)


def _ring_sum(values, r, edges):
    """Sum of real values and pixel counts in annuli defined by edges (DC pixel excluded)."""
    n = len(edges) - 1
    rr = r.ravel()
    idx = np.digitize(rr, edges) - 1
    ok = (idx >= 0) & (idx < n) & (rr > 0)
    s = np.bincount(idx[ok], weights=values.ravel()[ok], minlength=n)
    c = np.bincount(idx[ok], minlength=n)
    return s, c


def fourier_compare(a, b, window=True, nbins=64, nan_mode="interp", taper=16.0, interp_sigma=1.0):
    """2D power spectra, their log ratio, radial profiles and Fourier Ring Correlation.

    A, B and A-B share the same weight map: common validity mask (after optional
    inpainting), apodized by `taper` pixels, times the optional Hann window.
    """
    if a.shape != b.shape:
        raise ValueError(f"Shape mismatch: {a.shape} vs {b.shape}")
    if nan_mode not in NAN_MODES:
        raise ValueError(f"Unknown NaN mode '{nan_mode}'")
    if nan_mode == "interp":
        a, b = inpaint_common(a, b, interp_sigma)
    mask = np.isfinite(a) & np.isfinite(b)
    w = fft_weights(mask, window, taper)
    if not np.any(w > 0):
        raise ValueError("All FFT weights are zero (no valid pixel or taper too large)")
    norm = np.sum(w ** 2)  # Parseval: mean of P over the plane = weighted variance
    fa, fb, fd = (_centered_fft(x, w) for x in (a, b, a - b))
    pa, pb, pd = (np.abs(f) ** 2 / norm for f in (fa, fb, fd))
    fx, fy, r = _freq_grid(a.shape)
    edges = np.linspace(0.0, 0.5, nbins + 1)  # corners beyond Nyquist radius are excluded

    sa, cnt = _ring_sum(pa, r, edges)
    sb, _ = _ring_sum(pb, r, edges)
    sd, _ = _ring_sum(pd, r, edges)
    sx, _ = _ring_sum(np.real(fa * np.conj(fb)) / norm, r, edges)
    with np.errstate(divide="ignore", invalid="ignore"):
        prof = [np.where(cnt > 0, s / cnt, np.nan) for s in (sa, sb, sd)]
        frc = np.where(cnt > 0, sx / np.sqrt(sa * sb), np.nan)  # Fourier Ring Correlation

    return {
        "fx": fx, "fy": fy,
        "logps_a": _log10(pa), "logps_b": _log10(pb),
        "log_ratio": _log10(pa / np.where(pb > 0, pb, np.nan)),
        "f": 0.5 * (edges[1:] + edges[:-1]),
        "prof_a": prof[0], "prof_b": prof[1], "prof_d": prof[2],
        "frc": frc,
        "excluded_frac": 1.0 - mask.mean(),  # pixels without data in the FFT
    }
