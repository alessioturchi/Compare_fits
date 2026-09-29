"""Processing pipeline without GUI: settings, single-image view, comparison and batch comparison.

The GUI collects the settings, runs these functions in a worker thread and displays the
results; the same functions can be used from scripts.
"""
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import dataclass, fields

import numpy as np

from .fitscmp_analysis import difference, fourier_compare
from .fitscmp_calib import calibrate
from .fitscmp_io import load_image
from .fitscmp_registration import apply_shift, estimate_shift
from .fitscmp_tables import measure_spots

FRC_LEVEL = 0.5  # FRC level whose first crossing summarizes the Fourier comparison in the batch


class Cancelled(Exception):
    """The job was cancelled by the user."""


@dataclass
class Settings:
    """All processing parameters of one run (snapshot of the GUI options)."""
    normalization: str = "none"
    window: bool = True
    radial_bins: int = 64
    nan_mode: str = "interp"
    taper: float = 16.0
    interp_sigma: float = 1.0
    spot_method: str = "iso"
    saturation: float = None
    dark: np.ndarray = None
    dark_label: str = ""
    hot_filter: bool = False
    hot_nsigma: float = 5.0
    auto_shift: bool = False
    dy: float = 0.0
    dx: float = 0.0
    upsample: int = 100

    def describe(self):
        """JSON-friendly description (the dark is identified by its label)."""
        d = {f.name: getattr(self, f.name) for f in fields(self) if f.name != "dark"}
        d["dark_subtracted"] = self.dark is not None
        return d

    def calibration_text(self):
        parts = (["master dark subtracted"] if self.dark is not None else []) + (
            ["hot-pixel filter"] if self.hot_filter else [])
        return ", ".join(parts) or "none"


def load(refs, s):
    """Load and calibrate images: (images, saturation masks or None, hot-pixel counts)."""
    out = [calibrate(load_image(r), s.dark, s.hot_filter, s.hot_nsigma, s.saturation) for r in refs]
    return [o[0] for o in out], [o[1] for o in out], [o[2] for o in out]


def compare_images(a, b, s):
    """Registration (manual or automatic), difference and Fourier comparison of calibrated images."""
    if s.auto_shift:
        dy, dx = (round(v, 3) for v in estimate_shift(a, b, s.upsample, s.taper))  # as shown in the GUI
    else:
        dy, dx = float(s.dy), float(s.dx)
    b_shifted = apply_shift(b, dy, dx)
    a_n, b_n, d, stats = difference(a, b_shifted, s.normalization)
    fc = fourier_compare(a_n, b_n, s.window, s.radial_bins, s.nan_mode, s.taper, s.interp_sigma)
    return {"dy": dy, "dx": dx, "b_shifted": b_shifted, "d": d, "stats": stats, "fc": fc}


def run_single(refs, s):
    """Single-image view: calibrated images and their spot metrics."""
    imgs, sats, n_hot = load(refs, s)
    return {"imgs": imgs, "sats": sats, "n_hot": n_hot, "spots": measure_spots(imgs, s.spot_method, sats)}


def run_comparison(pair, s):
    """Comparison of A and B; spot metrics refer to the images as recorded (B before the shift)."""
    (a, b), sats, n_hot = load(pair, s)
    out = compare_images(a, b, s)
    out.update(a=a, b=b, sats=sats, n_hot=n_hot, spots=measure_spots([a, b], s.spot_method, sats))
    return out


def frc_crossing(f, frc, level=FRC_LEVEL):
    """First spatial frequency where the FRC drops below `level` (linear interpolation); NaN if never."""
    ok = np.isfinite(frc)
    f, frc = f[ok], frc[ok]
    below = np.nonzero(frc < level)[0]
    if below.size == 0 or below[0] == 0:
        return np.nan if below.size == 0 else float(f[0])
    k = below[0]
    return float(f[k - 1] + (level - frc[k - 1]) * (f[k] - f[k - 1]) / (frc[k] - frc[k - 1]))


# ---------------------------------------------------------------- batch
_SPOT_KEYS = ("background", "centroid_x", "centroid_y", "center_x", "center_y", "center_x_err", "center_y_err",
              "peak_value", "peak_x", "peak_y", "total_flux", "d_sigma_x", "d_sigma_y", "d_sigma",
              "fwhm_x", "fwhm_y", "phi_deg", "ellipticity", "n_saturated", "residual_ratio")


def _position(m):
    """Spot position of a metrics dict: fit center for the Gaussian fit, centroid otherwise."""
    return (m["center_x"], m["center_y"]) if m["method"] == "gauss" else (m["centroid_x"], m["centroid_y"])


def _spot_fields(m, err):
    row = {k: (float(m[k]) if m and k in m else np.nan) for k in _SPOT_KEYS}
    if m and m.get("round"):
        row["phi_deg"] = np.nan  # azimuth undefined for round beams (as in the side table)
    row["spot_error"] = err or ""
    return row


def batch_row(k, ref, a, a_spot, s):
    """Row of the batch table for image `ref` (index k) compared with the reference image `a`."""
    (b,), sats, n_hot = load([ref], s)
    m, err = measure_spots([b], s.spot_method, sats)[0]
    row = {"index": k, "image": ref.label(), "reference": False, "n_hot": n_hot[0], "error": ""}
    row.update(_spot_fields(m, err))
    ma = a_spot[0]
    if m and ma:  # displacement and flux ratio as recorded (before registration)
        (xa, ya), (xb, yb) = _position(ma), _position(m)
        row.update(dpos_x=xb - xa, dpos_y=yb - ya, flux_ratio=m["total_flux"] / ma["total_flux"])
    else:
        row.update(dpos_x=np.nan, dpos_y=np.nan, flux_ratio=np.nan)
    row.update(shift_dy=np.nan, shift_dx=np.nan, rms_rel=np.nan, pearson=np.nan, masked=np.nan, frc_f50=np.nan)
    if b.shape != a.shape:
        row["error"] = f"shape {b.shape} differs from the reference {a.shape}"
        return row
    c = compare_images(a, b, s)
    row.update(shift_dy=c["dy"], shift_dx=c["dx"], rms_rel=float(c["stats"]["rms_rel"]),
               pearson=float(c["stats"]["pearson"]), masked=float(c["stats"]["masked"]),
               frc_f50=frc_crossing(c["fc"]["f"], c["fc"]["frc"]))
    return row


_WORKER = {}  # per-process data of the batch workers: reference image, its spot, settings


def _worker_init(a, a_spot, s):
    _WORKER.update(a=a, a_spot=a_spot, s=s)


def _worker_task(k, ref):
    try:
        return batch_row(k, ref, _WORKER["a"], _WORKER["a_spot"], _WORKER["s"])
    except Exception as e:  # one failing image must not stop the batch
        return {"index": k, "image": ref.label(), "reference": False, "error": str(e)}


def run_batch(reference, refs, s, workers=1, progress=None, cancel=None):
    """Compare every image of `refs` with `reference`; rows sorted by index (0 = reference).

    workers > 1 uses a process pool (spawn start method: safe with a Tk application);
    progress(done, total) is called after each image; cancel() returning True stops the batch
    (raises Cancelled).
    """
    (a,), sats, n_hot = load([reference], s)
    a_spot = measure_spots([a], s.spot_method, sats)[0]
    ref_row = {"index": 0, "image": reference.label(), "reference": True, "n_hot": n_hot[0], "error": ""}
    ref_row.update(_spot_fields(*a_spot), dpos_x=0.0, dpos_y=0.0, flux_ratio=1.0 if a_spot[0] else np.nan,
                   shift_dy=0.0, shift_dx=0.0, rms_rel=0.0, pearson=1.0, masked=np.nan, frc_f50=np.nan)
    jobs = [(k, r) for k, r in enumerate(refs, start=1)]
    rows, done = [ref_row], 0
    if workers <= 1:
        _worker_init(a, a_spot, s)
        for k, r in jobs:
            if cancel and cancel():
                raise Cancelled()
            rows.append(_worker_task(k, r))
            done += 1
            if progress:
                progress(done, len(jobs))
    else:
        with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn"),
                                 initializer=_worker_init, initargs=(a, a_spot, s)) as pool:
            futures = [pool.submit(_worker_task, k, r) for k, r in jobs]
            for fut in as_completed(futures):
                if cancel and cancel():
                    pool.shutdown(wait=False, cancel_futures=True)
                    raise Cancelled()
                rows.append(fut.result())
                done += 1
                if progress:
                    progress(done, len(jobs))
    return sorted(rows, key=lambda r: r["index"])
