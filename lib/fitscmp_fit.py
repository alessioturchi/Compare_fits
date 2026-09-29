"""Elliptical 2D Gaussian fit of a laser spot, initialized by the ISO 11146 analysis."""
import numpy as np
from scipy.optimize import least_squares

from .fitscmp_spot import ROUND_LIMIT, iso_spot, saturated_mask

FWHM_PER_SIGMA = 2.0 * np.sqrt(2.0 * np.log(2.0))


def _model(p, x, y):
    """Constant background + elliptical Gaussian; s1 is the sigma along the axis at angle theta."""
    bg, amp, x0, y0, s1, s2, th = p
    c, s = np.cos(th), np.sin(th)
    u, v = (x - x0) * c + (y - y0) * s, -(x - x0) * s + (y - y0) * c
    return bg + amp * np.exp(-0.5 * ((u / s1) ** 2 + (v / s2) ** 2))


def _iso_axes(s1, s2, th, e1, e2):
    """Map (s1, s2, theta) to the ISO convention: axis closer to x first, |phi| <= 45 deg."""
    th = (th + np.pi / 2) % np.pi - np.pi / 2  # (-90, 90] deg
    if abs(th) <= np.pi / 4:
        return s1, s2, th, e1, e2
    return s2, s1, th - np.sign(th) * np.pi / 2, e2, e1


def gauss_spot(img, saturation=None):
    """Least-squares fit on the bounding box of the ISO integration area.

    Saturated pixels (`saturation`: level or boolean mask, see saturated_mask) and NaN pixels
    are excluded from the fit.
    Uncertainties are 1-sigma, from the heteroscedasticity-consistent (sandwich) covariance
    (White 1980), which does not need a noise model.
    """
    iso = iso_spot(img, saturation)
    ny, nx = img.shape
    au, av, ap = iso["area"]
    hw = 1.5 * (abs(au * np.cos(ap)) + abs(av * np.sin(ap)))  # half-sizes of the area bounding box
    hh = 1.5 * (abs(au * np.sin(ap)) + abs(av * np.cos(ap)))
    x_lo, x_hi = max(int(iso["centroid_x"] - hw), 0), min(int(iso["centroid_x"] + hw) + 1, nx)
    y_lo, y_hi = max(int(iso["centroid_y"] - hh), 0), min(int(iso["centroid_y"] + hh) + 1, ny)
    cut = img[y_lo:y_hi, x_lo:x_hi]
    yy, xx = np.mgrid[y_lo:y_hi, x_lo:x_hi].astype(float)
    sat = saturated_mask(img, saturation)
    use = np.isfinite(cut) if sat is None else np.isfinite(cut) & ~sat[y_lo:y_hi, x_lo:x_hi]
    x, y, z = xx[use], yy[use], cut[use]

    s_x, s_y = iso["d_sigma_x"] / 4.0, iso["d_sigma_y"] / 4.0  # ISO widths as starting sigmas
    p0 = [iso["background"], max(np.max(z) - iso["background"], 1e-6), iso["centroid_x"], iso["centroid_y"],
          s_x, s_y, np.radians(iso["phi_deg"])]
    lo = [-np.inf, 0.0, x_lo - 0.5, y_lo - 0.5, 0.3, 0.3, -np.pi]
    hi = [np.inf, np.inf, x_hi - 0.5, y_hi - 0.5, np.inf, np.inf, np.pi]
    p0 = np.clip(p0, np.nextafter(lo, hi), np.nextafter(hi, lo))
    res = least_squares(lambda p: _model(p, x, y) - z, p0, bounds=(lo, hi), x_scale="jac")

    resid_var = np.sum(res.fun ** 2) / max(z.size - len(p0), 1)
    try:  # heteroscedasticity-consistent (sandwich) covariance: shot noise varies across the spot
        j = res.jac
        bread = np.linalg.inv(j.T @ j)
        meat = (j * (res.fun ** 2)[:, None]).T @ j
        dof_corr = z.size / max(z.size - len(p0), 1)
        err = np.sqrt(np.diag(bread @ meat @ bread) * dof_corr)
    except np.linalg.LinAlgError:
        err = np.full(len(p0), np.nan)
    bg, amp, x0, y0, sig1, sig2, th = res.x
    sx_, sy_, phi, ex, ey = _iso_axes(sig1, sig2, th, err[4], err[5])
    ell = min(sx_, sy_) / max(sx_, sy_)
    out = {"method": "gauss", "background": bg, "amplitude": amp, "center_x": x0, "center_y": y0,
           "center_x_err": err[2], "center_y_err": err[3],
           "d_sigma_x": 4 * sx_, "d_sigma_y": 4 * sy_, "d_sigma_x_err": 4 * ex, "d_sigma_y_err": 4 * ey,
           "d_sigma": 2 * np.sqrt(2) * np.hypot(sx_, sy_),
           "fwhm_x": FWHM_PER_SIGMA * sx_, "fwhm_y": FWHM_PER_SIGMA * sy_,
           "phi_deg": np.degrees(phi), "ellipticity": ell, "round": bool(ell > ROUND_LIMIT),
           "total_flux": 2 * np.pi * amp * sig1 * sig2, "residual_rms": np.sqrt(resid_var),
           "residual_ratio": np.sqrt(resid_var) / iso["background_rms"],  # >> 1: non-Gaussian profile
           "n_fit": int(z.size), "converged": bool(res.success),
           "fit_box": (x_lo, x_hi, y_lo, y_hi)}
    out.update({k: iso[k] for k in ("peak_value", "peak_x", "peak_y", "n_peak") + (
        ("n_saturated",) if saturation is not None else ())})
    return out
