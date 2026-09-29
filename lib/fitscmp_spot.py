"""Laser-spot metrics: ISO 11146 second-moment analysis and the Rainer (p10) method.

Coordinates are 0-based pixel indices (x = column, y = row) with integer values at
pixel centers, as in the image panels.
"""
import numpy as np
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

SPOT_METHODS = {"iso": "ISO 11146", "gauss": "Gaussian fit", "rainer": "Rainer (p10)"}
ROUND_LIMIT = 0.87        # ISO 11146: ellipticity above this value -> round beam
RAINER_FWHM = 2.3548      # FWHM / sigma factor used by the Rainer method (Gaussian profile)
_MIN_OUTSIDE = 100        # minimum number of pixels outside the integration area for the background
_DETECT_NSIGMA = 3.0      # detection threshold for the initial estimate [background rms]


class SpotError(ValueError):
    """The spot cannot be measured (no positive signal, degenerate moments)."""


def saturated_mask(img, saturation):
    """Saturated pixels: `saturation` is None, a level (pixels >= level) or a boolean mask.

    A mask lets the caller evaluate saturation on the raw values, e.g. before dark subtraction.
    """
    if saturation is None:
        return None
    if isinstance(saturation, np.ndarray) and saturation.dtype == bool:
        return saturation
    return np.isfinite(img) & (img >= saturation)


def _peak(img, finite, saturation):
    """Peak value and position (mean position of all pixels at the maximum) and saturation count."""
    vmax = np.max(img[finite])
    ys, xs = np.nonzero(finite & (img == vmax))
    out = {"peak_value": vmax, "peak_x": xs.mean(), "peak_y": ys.mean(), "n_peak": xs.size}
    sat = saturated_mask(img, saturation)
    if sat is not None:
        out["n_saturated"] = int(np.sum(finite & sat))
    return out


def _moments(e, xx, yy):
    """Total, first moments and second central moments of the weight map e."""
    t = e.sum()
    if not t > 0:
        raise SpotError("No positive signal above the background")
    xc, yc = (e * xx).sum() / t, (e * yy).sum() / t
    dx, dy = xx - xc, yy - yc
    return t, xc, yc, (e * dx * dx).sum() / t, (e * dy * dy).sum() / t, (e * dx * dy).sum() / t


def iso_widths(sxx, syy, sxy):
    """ISO 11146-1 beam widths d_sigma_x, d_sigma_y along the principal axes and azimuth phi [rad].

    phi is the angle between the x axis and the principal axis closer to it (|phi| <= 45 deg);
    d_sigma_x is the width along that axis, d_sigma_y the width along the other one.
    """
    d = sxx - syy
    root = np.sqrt(d * d + 4.0 * sxy * sxy)
    g = np.sign(d) if d != 0 else 1.0
    dsx = 2.0 * np.sqrt(2.0) * np.sqrt(sxx + syy + g * root)
    dsy = 2.0 * np.sqrt(2.0) * np.sqrt(max(sxx + syy - g * root, 0.0))
    if d != 0:
        phi = 0.5 * np.arctan(2.0 * sxy / d)
    else:
        phi = np.sign(sxy) * np.pi / 4.0  # 0 for a round beam
    return dsx, dsy, phi


def _detect(img, finite, bg, noise):
    """Initial spot region: connected component with the largest flux among the pixels whose
    3x3 median-filtered value exceeds bg + 3 rms (isolated hot pixels do not pass the filter)."""
    smooth = ndimage.median_filter(np.where(finite, img, bg), size=3)
    labels, n = ndimage.label(finite & (smooth > bg + _DETECT_NSIGMA * noise))
    if n == 0:
        raise SpotError(f"No pixel above background + {_DETECT_NSIGMA:g} rms")
    flux = ndimage.sum(np.where(finite, img - bg, 0.0), labels, index=np.arange(1, n + 1))
    return labels == (np.argmax(flux) + 1)


def _area(xx, yy, xc, yc, dsx, dsy, phi):
    """ISO integration area: rectangle of sides 3 d_sigma_x, 3 d_sigma_y aligned with the principal axes."""
    u = (xx - xc) * np.cos(phi) + (yy - yc) * np.sin(phi)
    v = -(xx - xc) * np.sin(phi) + (yy - yc) * np.cos(phi)
    return (np.abs(u) <= 1.5 * dsx) & (np.abs(v) <= 1.5 * dsy)


def iso_spot(img, saturation=None, max_iter=50):
    """ISO 11146 spot analysis with iterative background and integration area.

    Start: sigma-clipped median of the whole frame as background, and moments of the detected
    spot region (see _detect). Then background = sigma-clipped mean of the pixels outside the
    integration area. Integration area: rectangle centered on
    the centroid, aligned with the principal axes, with sides 3 x d_sigma_x and 3 x d_sigma_y;
    for round beams (azimuth undefined) a square of side 3 x d_sigma aligned with the frame.
    Iterated until the set of pixels in the area is stationary or repeats (discrete limit cycle).
    """
    finite = np.isfinite(img)
    if not finite.any():
        raise SpotError("No finite pixel")
    ny, nx = img.shape
    yy, xx = np.mgrid[:ny, :nx].astype(float)
    _, bg, noise = sigma_clipped_stats(img[finite], sigma=3.0, maxiters=10)  # robust start
    inside, seen, converged = _detect(img, finite, bg, noise), set(), False
    for it in range(1, max_iter + 1):
        t, xc, yc, sxx, syy, sxy = _moments(np.where(inside, img - bg, 0.0), xx, yy)
        if sxx <= 0 or syy <= 0:
            raise SpotError("degenerate second moments: no isolated spot, or spot larger than the frame")
        dsx, dsy, phi = iso_widths(sxx, syy, sxy)
        ds = 2.0 * np.sqrt(2.0) * np.sqrt(sxx + syy)
        is_round = min(dsx, dsy) / max(dsx, dsy) > ROUND_LIMIT
        area = (ds, ds, 0.0) if is_round else (dsx, dsy, phi)
        new = finite & _area(xx, yy, xc, yc, *area)
        key = np.packbits(new).tobytes()
        if key in seen:  # stationary area or limit cycle: further iterations change nothing
            converged = True
            break
        seen.add(key)
        inside = new
        outside = finite & ~inside
        if outside.sum() >= _MIN_OUTSIDE:
            bg, _, noise = sigma_clipped_stats(img[outside], sigma=3.0, maxiters=10, cenfunc="mean")

    # corners of the final integration area, to check that it lies within the frame
    au, av, ap = area
    cu, cv = np.array([1, 1, -1, -1]) * 1.5 * au, np.array([1, -1, 1, -1]) * 1.5 * av
    cx = xc + cu * np.cos(ap) - cv * np.sin(ap)
    cy = yc + cu * np.sin(ap) + cv * np.cos(ap)
    in_frame = bool(np.all((cx >= -0.5) & (cx <= nx - 0.5) & (cy >= -0.5) & (cy <= ny - 0.5)))
    out = {"method": "iso", "background": bg, "background_rms": noise,
           "centroid_x": xc, "centroid_y": yc, "total_flux": t,
           "d_sigma_x": dsx, "d_sigma_y": dsy, "d_sigma": ds, "phi_deg": np.degrees(phi),
           "ellipticity": min(dsx, dsy) / max(dsx, dsy), "n_iter": it, "converged": converged,
           "round": bool(is_round), "area": area, "area_in_frame": in_frame,
           "n_nan_in_area": int(np.sum(~finite & _area(xx, yy, xc, yc, *area)))}
    out.update(_peak(img, finite, saturation))
    return out


def rainer_spot(img, saturation=None):
    """Rainer (p10) method, as in the Acquisition Manager camera view.

    Background = 10th percentile of the frame; negative residuals clipped to zero; centroid
    and sigma from the x and y projections of the whole frame; FWHM = 2.3548 sigma.
    """
    finite = np.isfinite(img)
    if not finite.any():
        raise SpotError("No finite pixel")
    bg = float(np.percentile(img[finite], 10))
    e = np.clip(np.where(finite, img - bg, 0.0), 0.0, None)
    total = e.sum()
    if total < 1.0:  # same threshold as the original implementation
        raise SpotError("No signal above the 10th-percentile background")
    x, y = np.arange(img.shape[1], dtype=float), np.arange(img.shape[0], dtype=float)
    px, py = e.sum(axis=0), e.sum(axis=1)
    xc, yc = (px * x).sum() / total, (py * y).sum() / total
    sx = np.sqrt((px * (x - xc) ** 2).sum() / total)
    sy = np.sqrt((py * (y - yc) ** 2).sum() / total)
    out = {"method": "rainer", "background": bg, "centroid_x": xc, "centroid_y": yc,
           "total_flux": total, "fwhm_x": RAINER_FWHM * sx, "fwhm_y": RAINER_FWHM * sy}
    out.update(_peak(img, finite, saturation))
    return out


def spot_metrics(img, method="iso", saturation=None):
    """Dispatch to the selected method ('iso', 'gauss' or 'rainer')."""
    if method == "iso":
        return iso_spot(img, saturation)
    if method == "gauss":
        from .fitscmp_fit import gauss_spot  # local import: fitscmp_fit depends on this module
        return gauss_spot(img, saturation)
    if method == "rainer":
        return rainer_spot(img, saturation)
    raise ValueError(f"Unknown spot method '{method}'")
