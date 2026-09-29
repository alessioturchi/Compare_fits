"""Validation of lib/fitscmp_pipeline.py (batch, process pool, cancel) and lib/fitscmp_export.py.

Run with `python -m pytest tests` or `python tests/test_pipeline.py` from the project root.
"""
import importlib.util
import json
import os
import sys
import tempfile

import numpy as np
from astropy.io import fits

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lib.fitscmp_export import export_batch, export_result, summary_rows  # noqa: E402
from lib.fitscmp_io import list_images  # noqa: E402
from lib.fitscmp_pipeline import Cancelled, Settings, frc_crossing, run_batch, run_comparison  # noqa: E402
from lib.fitscmp_tables import batch_columns, batch_values, build_table  # noqa: E402

TMP = tempfile.mkdtemp(prefix="fitscmp_test_")
N = 8


def drift_file():
    """Cube of N frames: spot moving by (+0.5, -0.25) pix per frame, width growing 1% per frame."""
    path = os.path.join(TMP, "drift.fits")
    if not os.path.exists(path):
        rng = np.random.default_rng(3)
        yy, xx = np.mgrid[:200, :240].astype(float)
        planes = []
        for k in range(N):
            s = 12 * (1 + 0.01 * k)
            sig = 2000 * np.exp(-((xx - 110 - 0.5 * k) ** 2 + (yy - 100 + 0.25 * k) ** 2) / (2 * s ** 2))
            planes.append(sig + 64 + rng.normal(0, 5, sig.shape) + rng.normal(0, np.sqrt(sig / 2)))
        fits.PrimaryHDU(np.array(planes, np.float32)).writeto(path)
    return list_images(path)


def test_frc_crossing():
    f = np.array([0.1, 0.2, 0.3, 0.4])
    assert np.isclose(frc_crossing(f, np.array([1, 0.8, 0.4, 0.1])), 0.275)  # linear interpolation
    assert np.isnan(frc_crossing(f, np.array([1, 0.9, 0.8, 0.7])))
    assert frc_crossing(f, np.array([0.3, 0.2, 0.1, 0.0])) == 0.1


def test_batch_recovers_drift():
    refs = drift_file()
    rows = run_batch(refs[0], refs[1:], Settings())
    k = np.arange(N)
    assert [r["index"] for r in rows] == list(k) and rows[0]["reference"]
    assert np.max(np.abs(np.array([r["dpos_x"] for r in rows]) - 0.5 * k)) < 0.1
    assert np.max(np.abs(np.array([r["dpos_y"] for r in rows]) + 0.25 * k)) < 0.1
    assert np.max(np.abs(np.array([r["d_sigma"] for r in rows]) / (48 * (1 + 0.01 * k)) - 1)) < 0.01
    assert np.allclose([r["flux_ratio"] for r in rows], (1 + 0.01 * k) ** 2, rtol=0.01)  # fixed amplitude
    assert all(np.isnan(r["phi_deg"]) for r in rows)  # round spots: azimuth undefined


def test_process_pool_matches_sequential():
    refs = drift_file()
    s = Settings(auto_shift=True)
    seq, par = run_batch(refs[0], refs[1:], s, workers=1), run_batch(refs[0], refs[1:], s, workers=2)
    for a, b in zip(seq, par):
        assert a.keys() == b.keys()
        for key in a:
            assert (a[key] == b[key]) or (isinstance(a[key], float) and np.isnan(a[key]) and np.isnan(b[key]))


def test_cancel():
    refs = drift_file()
    for workers in (1, 2):
        done = []
        try:
            run_batch(refs[0], refs[1:], Settings(), workers, progress=lambda d, n: done.append(d),
                      cancel=lambda: len(done) >= 2)
        except Cancelled:
            assert len(done) < N - 1
            continue
        raise AssertionError("batch not cancelled")


def test_exports_round_trip():
    refs = drift_file()
    s = Settings(saturation=4095)
    res = run_comparison((refs[0], refs[3]), s)
    rows, note = build_table([res["a"], res["b"]], ["A", "B"], res["spots"], "iso", 4095, delta=True)
    table = (["A", "B", "B \u2212 A"], rows, note)
    fmts = ["csv"] + (["parquet"] if all(importlib.util.find_spec(m) for m in ("pandas", "pyarrow")) else [])
    for fmt in fmts:
        paths = export_result(os.path.join(TMP, f"cmp_{fmt}"), table, {"settings": s.describe(), "stats": res["stats"]},
                              res["fc"], fmt)
        assert len(paths) == 3 and all(os.path.exists(p) for p in paths)
        with open(paths[2]) as f:
            meta = json.load(f)
        assert meta["settings"]["saturation"] == 4095 and meta["files"][0].startswith(f"cmp_{fmt}_summary")
    tidy = summary_rows(table)
    dx = [r for r in tidy if r[1] == "Centroid x" and r[3] == "B \u2212 A"][0][4]
    assert np.isclose(dx, res["spots"][1][0]["centroid_x"] - res["spots"][0][0]["centroid_x"])
    ratio = [r for r in tidy if r[1] == "Total flux" and r[3] == "B \u2212 A"][0]
    assert ratio[5] == "ratio B/A" and abs(ratio[4] - 1.0609) < 0.01
    brows = run_batch(refs[0], refs[1:3], s)
    cols = batch_columns("iso", 4095)
    paths = export_batch(os.path.join(TMP, "batch"), [c[2] for c in cols], batch_values(brows, cols), {"n": 3})
    with open(paths[0]) as f:
        header = f.readline().strip().split(",")
    assert header[:4] == ["index", "image", "centroid_x_pix", "centroid_y_pix"] and "frc_f50_cycles_per_pix" in header


if __name__ == "__main__":
    tests = [(k, v) for k, v in dict(globals()).items() if k.startswith("test_")]
    for name, fn in tests:
        fn()
        print(f"ok  {name}")
    print(f"{len(tests)} tests passed")
