"""Validation of lib/fitscmp_calib.py (master dark, hot pixels) on simulated camera frames.

Run with `python -m pytest tests` or `python tests/test_calib.py` from the project root.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.fitscmp_calib import calibrate, hot_pixels  # noqa: E402
from lib.fitscmp_spot import iso_spot  # noqa: E402

NY, NX = 512, 640
YY, XX = np.mgrid[:NY, :NX].astype(float)


def spot(x0=300.0, y0=260.0, s=25.0, peak=2000.0):
    return peak * np.exp(-((XX - x0) ** 2 + (YY - y0) ** 2) / (2 * s ** 2))


def noisy(signal, rng, offset=64.0, ron=5.0):
    return signal + offset + rng.normal(0, ron, signal.shape) + rng.normal(0, np.sqrt(np.maximum(signal, 0) / 2))


def add_hot(img, rng, n=30, lo=50.0, hi=4000.0):
    """Isolated hot pixels at random positions (at least 2 pixels apart from each other)."""
    out, pos = img.copy(), set()
    while len(pos) < n:
        y, x = int(rng.integers(1, NY - 1)), int(rng.integers(1, NX - 1))
        if all(abs(y - py) > 2 or abs(x - px) > 2 for py, px in pos):
            pos.add((y, x))
    ys, xs = np.array(sorted(pos)).T
    out[ys, xs] += rng.uniform(lo, hi, n)
    return out, (ys, xs)


def test_no_false_positives_on_noise_and_spot():
    rng = np.random.default_rng(0)
    big = np.random.default_rng(1).normal(64, 5, (1024, 1280))
    assert hot_pixels(big).sum() == 0                     # 1.3 Mpix of pure noise
    assert hot_pixels(noisy(spot(), rng)).sum() == 0      # far-field spot


def test_sharp_spot_peak_is_not_flagged():
    rng = np.random.default_rng(2)
    for s in (1.0, 1.5, 3.0):                            # down to FWHM = 2.4 pix
        assert hot_pixels(noisy(spot(300.3, 260.6, s, 4000.0), rng)).sum() == 0


def test_hot_pixels_are_found():
    rng = np.random.default_rng(3)
    img, (ys, xs) = add_hot(noisy(spot(), rng), rng, lo=50.0)  # down to 10 x read noise
    mask = hot_pixels(img)
    assert mask[ys, xs].all() and mask.sum() == len(ys)


def test_filter_removes_the_width_bias():
    rng = np.random.default_rng(4)
    clean = noisy(spot(), rng)
    hot, _ = add_hot(clean, rng)
    ref = iso_spot(clean)["d_sigma"]
    assert iso_spot(hot)["d_sigma"] / ref - 1 > 0.01      # biased without the filter
    fixed, _, n_hot = calibrate(hot, hot_filter=True)
    assert n_hot == 30 and abs(iso_spot(fixed)["d_sigma"] / ref - 1) < 0.002


def test_master_dark_and_raw_saturation():
    rng = np.random.default_rng(5)
    pattern = 64.0 + 3.0 * np.sin(XX / 40.0)             # offset with fixed pattern
    pattern, _ = add_hot(pattern, rng, lo=200.0, hi=1500.0)  # warm pixels, present in every frame
    dark = np.mean([pattern + rng.normal(0, 5, pattern.shape) for _ in range(16)], axis=0)
    raw = np.minimum(noisy(spot(peak=8000.0), rng, offset=0.0) + pattern, 4095.0)  # 12-bit saturation
    cal, sat, _ = calibrate(raw, dark=dark, saturation=4095)
    r = iso_spot(cal, saturation=sat)
    assert abs(r["background"]) < 0.1                     # offset and pattern removed
    assert r["n_saturated"] == int((raw >= 4095).sum()) > 0
    assert iso_spot(cal, saturation=4095).get("n_saturated") == 0  # level on dark-subtracted data misses them
    ok = iso_spot(np.minimum(noisy(spot(peak=2000.0), rng, offset=0.0) + pattern, 4095.0) - dark)
    assert abs(ok["d_sigma"] / 100 - 1) < 0.005           # warm pixels removed by the dark


def test_dark_shape_mismatch_raises():
    try:
        calibrate(np.zeros((10, 10)), dark=np.zeros((10, 11)))
    except ValueError:
        return
    raise AssertionError("shape mismatch not detected")


if __name__ == "__main__":
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        fn()
        print(f"ok  {name}")
    print(f"{len(tests)} tests passed")
