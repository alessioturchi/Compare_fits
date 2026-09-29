"""Validation of lib/fitscmp_spot.py on simulated laser spots.

Run with `python -m pytest tests` or `python tests/test_spot.py` from the project root.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.fitscmp_fit import gauss_spot  # noqa: E402
from lib.fitscmp_spot import SpotError, iso_spot, rainer_spot  # noqa: E402

NY, NX = 512, 640
OFFSET, RON = 64.0, 5.0  # camera dark offset and read noise [ADU]
YY, XX = np.mgrid[:NY, :NX].astype(float)


def camera(signal, seed=0):
    """Signal + dark offset + read noise + shot noise (gain 2 e-/ADU)."""
    rng = np.random.default_rng(seed)
    return signal + OFFSET + rng.normal(0, RON, signal.shape) + rng.normal(0, np.sqrt(np.maximum(signal, 0) / 2))


def gaussian(x0, y0, s_major, s_minor, theta_deg, peak=2000.0):
    """Elliptical Gaussian with the major axis at theta from the x axis."""
    t = np.radians(theta_deg)
    u = (XX - x0) * np.cos(t) + (YY - y0) * np.sin(t)
    v = -(XX - x0) * np.sin(t) + (YY - y0) * np.cos(t)
    return peak * np.exp(-0.5 * ((u / s_major) ** 2 + (v / s_minor) ** 2))


def rainer_reference(gray):
    """Verbatim copy of CameraView.compute_centroid (Acquisition Manager, lib/camera_view_lib.py)."""
    bg = float(np.percentile(gray, 10))
    img = gray.astype(np.float64) - bg
    np.clip(img, 0, None, out=img)
    total = img.sum()
    if total < 1.0:
        return None
    h, w = img.shape
    xx = np.arange(w, dtype=np.float64)
    yy = np.arange(h, dtype=np.float64)
    px = img.sum(axis=0)
    py = img.sum(axis=1)
    xc = (px * xx).sum() / total
    yc = (py * yy).sum() / total
    sx = np.sqrt(((px * (xx - xc) ** 2).sum()) / total)
    sy = np.sqrt(((py * (yy - yc) ** 2).sum()) / total)
    return xc, yc, 2.3548 * sx, 2.3548 * sy


def test_round_gaussian():
    r = iso_spot(camera(gaussian(250.3, 200.7, 20, 20, 0)))
    assert abs(r["centroid_x"] - 250.3) < 0.05 and abs(r["centroid_y"] - 200.7) < 0.05
    assert abs(r["d_sigma"] / 80.0 - 1) < 0.005          # d_sigma = 4 sigma
    assert r["round"] and r["converged"] and r["area_in_frame"]
    assert abs(r["background"] - OFFSET) < 0.05
    assert abs(r["total_flux"] / (2000 * 2 * np.pi * 20 * 20) - 1) < 0.005


def test_elliptical_major_closer_to_x():
    r = iso_spot(camera(gaussian(320, 256, 30, 15, 30)))
    assert abs(r["d_sigma_x"] / 120 - 1) < 0.01 and abs(r["d_sigma_y"] / 60 - 1) < 0.01
    assert abs(r["phi_deg"] - 30) < 0.5 and abs(r["ellipticity"] - 0.5) < 0.01 and not r["round"]


def test_elliptical_major_closer_to_y():
    # major axis at 60 deg: the principal axis closer to x is the minor one, at -30 deg
    r = iso_spot(camera(gaussian(320, 256, 30, 15, 60)))
    assert abs(r["phi_deg"] + 30) < 0.5
    assert abs(r["d_sigma_x"] / 60 - 1) < 0.01 and abs(r["d_sigma_y"] / 120 - 1) < 0.01


def test_top_hat_near_field():
    # uniform disk of radius R: sigma_x = sigma_y = R/2, so d_sigma = 2R (the diameter)
    disk = 1500.0 * (((XX - 330) ** 2 + (YY - 250) ** 2) < 100 ** 2)
    r = iso_spot(camera(disk))
    assert abs(r["d_sigma"] / 200 - 1) < 0.01 and r["round"]
    assert abs(r["centroid_x"] - 330) < 0.05 and abs(r["centroid_y"] - 250) < 0.05


def test_hot_pixels_outside_area_do_not_bias():
    # hot pixels inside the integration area DO bias the second moments (dark subtraction needed);
    # outside it they must not affect the sigma-clipped background
    frame = camera(gaussian(300, 260, 25, 25, 0))
    ref = iso_spot(frame)  # integration area: x in [150, 450], y in [110, 410]
    hot = frame.copy()
    rng = np.random.default_rng(5)
    xs = np.concatenate([rng.integers(0, 140, 15), rng.integers(460, NX, 15)])
    hot[rng.integers(0, NY, 30), xs] = 4095.0  # 12-bit saturated hot pixels
    r = iso_spot(hot)
    assert abs(r["background"] - ref["background"]) < 0.05
    assert abs(r["centroid_x"] - ref["centroid_x"]) < 0.1 and abs(r["d_sigma"] / ref["d_sigma"] - 1) < 0.01


def test_area_outside_frame_is_flagged():
    r = iso_spot(camera(gaussian(40, 256, 25, 25, 0)))  # 3 x d_sigma = 300 pix does not fit
    assert not r["area_in_frame"]


def test_saturated_plateau_peak_position():
    frame = np.minimum(camera(gaussian(300.0, 200.0, 25, 25, 0, peak=8000.0)), 4095.0)
    r = iso_spot(frame, saturation=4095)
    assert r["n_saturated"] == r["n_peak"] > 1
    assert abs(r["peak_x"] - 300) < 0.5 and abs(r["peak_y"] - 200) < 0.5


def test_nan_pixels_are_ignored():
    frame = camera(gaussian(300, 260, 25, 25, 0))
    frame[10:20, 10:600] = np.nan  # dead rows far from the spot
    r = iso_spot(frame)
    assert abs(r["centroid_x"] - 300) < 0.05 and abs(r["d_sigma"] / 100 - 1) < 0.005


def test_constant_frame_raises():
    for f in (iso_spot, rainer_spot):
        try:
            f(np.full((64, 64), 100.0))
        except SpotError:
            continue
        raise AssertionError(f"{f.__name__} did not raise")


def test_rainer_matches_acquisition_manager():
    gray = np.round(camera(gaussian(250, 180, 20, 12, 25))).astype(np.uint16)
    xc, yc, fx, fy = rainer_reference(gray)
    r = rainer_spot(gray.astype(float))
    assert np.allclose([r["centroid_x"], r["centroid_y"], r["fwhm_x"], r["fwhm_y"]], [xc, yc, fx, fy],
                       rtol=1e-12, atol=0)


def test_fit_elliptical_both_conventions():
    for theta, (dx, dy, phi) in ((30, (120, 60, 30)), (60, (60, 120, -30))):
        r = gauss_spot(camera(gaussian(320, 256, 30, 15, theta)))
        assert abs(r["d_sigma_x"] / dx - 1) < 0.002 and abs(r["d_sigma_y"] / dy - 1) < 0.002
        assert abs(r["phi_deg"] - phi) < 0.2 and abs(r["center_x"] - 320) < 0.02
        assert abs(r["total_flux"] / (2000 * 2 * np.pi * 30 * 15) - 1) < 0.002 and r["converged"]


def test_fit_excludes_saturated_pixels():
    frame = np.minimum(camera(gaussian(300.0, 200.0, 25, 25, 0, peak=8000.0)), 4095.0)
    r = gauss_spot(frame, saturation=4095)
    assert abs(r["d_sigma"] / 100 - 1) < 0.002 and abs(r["amplitude"] / 8000 - 1) < 0.01
    assert iso_spot(frame, saturation=4095)["d_sigma"] > 104  # moments are inflated by the plateau


def test_fit_is_robust_to_hot_pixels():
    frame = camera(gaussian(300, 260, 25, 25, 0))
    hot = frame.copy()
    rng = np.random.default_rng(5)
    hot[rng.integers(0, NY, 30), rng.integers(0, NX, 30)] = 4095.0
    assert abs(gauss_spot(hot)["d_sigma"] / gauss_spot(frame)["d_sigma"] - 1) < 0.002


def test_fit_uncertainties_are_calibrated():
    pulls = []
    for seed in range(20):
        r = gauss_spot(camera(gaussian(320.3, 256.6, 28, 22, 0, peak=200.0), seed))
        pulls.append([(r["center_x"] - 320.3) / r["center_x_err"], (r["d_sigma_x"] - 112) / r["d_sigma_x_err"]])
    assert np.all((np.std(pulls, axis=0) > 0.6) & (np.std(pulls, axis=0) < 1.5))


def test_fit_flags_non_gaussian_profile():
    disk = 1500.0 * (((XX - 330) ** 2 + (YY - 250) ** 2) < 100 ** 2)
    assert gauss_spot(camera(disk))["residual_ratio"] > 10
    assert gauss_spot(camera(gaussian(300, 260, 25, 25, 0)))["residual_ratio"] < 3


if __name__ == "__main__":
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        fn()
        print(f"ok  {name}")
    print(f"{len(tests)} tests passed")
