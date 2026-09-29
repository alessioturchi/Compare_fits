import os, sys, time, glob, warnings
sys.path.insert(0, os.getcwd())
import numpy as np
from astropy.io import fits

def make_drift(path, n=12):
    """n frames: spot drifting by (+0.5, -0.25) pix per frame, width growing 1%/frame."""
    rng = np.random.default_rng(3); yy, xx = np.mgrid[:240, :320].astype(float)
    planes = []
    for k in range(n):
        s = 12 * (1 + 0.01 * k); x0, y0 = 150 + 0.5 * k, 120 - 0.25 * k
        sig = 2000 * np.exp(-((xx - x0) ** 2 + (yy - y0) ** 2) / (2 * s ** 2))
        planes.append(sig + 64 + rng.normal(0, 5, sig.shape) + rng.normal(0, np.sqrt(sig / 2)))
    fits.PrimaryHDU(np.array(planes, np.float32)).writeto(path, overwrite=True)

if __name__ == "__main__":
    warnings.simplefilter("error")
    import fitscmp_gui as g
    make_drift("/tmp/drift.fits")
    infos = []; g.messagebox.showinfo = lambda *a, **k: infos.append(a[1])
    app = g.App(); app.withdraw(); app.paths = ["/tmp/drift.fits"]
    def wait():
        ticks = [0]
        def tick():
            if app._job is not None: ticks[0] += 1; app.after(20, tick)
        app.after(20, tick)
        while app._job is not None: app.update(); time.sleep(0.01)
        return ticks[0]
    last = lambda cls: [w for w in app.winfo_children() if isinstance(w, cls)][-1]
    def batch(ref_index=0):
        def f():
            d = [x for x in app.winfo_children() if isinstance(x, g.BatchDialog)][0]; d._cb.current(ref_index); d._ok()
        app.after(300, f); app.batch()
    # 1) batch, sequential and with 2 processes: identical, drift recovered
    res = {}
    for w in (1, 2):
        app.cfg["batch_workers"] = w; batch(); t = time.perf_counter(); ticks = wait(); bw = last(g.BatchWindow)
        res[w] = (bw.values, bw.columns); print(f"1) workers={w}: {time.perf_counter() - t:.1f} s, GUI ticks during job: {ticks}, status: '{app.status.get()}'")
    print("   identical tables:", all(str(a) == str(b) for a, b in zip(res[1][0], res[2][0])))
    cols = res[1][1]; vals = np.array([[v if isinstance(v, (int, float)) and not isinstance(v, bool) else np.nan for v in r] for r in res[1][0]], float)
    ix, iy, idn = cols.index("Δx B−A [pix]"), cols.index("Δy B−A [pix]"), cols.index("dσ [pix]")
    k = np.arange(12)
    print(f"   Δx - 0.5k max |err| {np.max(np.abs(vals[:, ix] - 0.5 * k)):.4f} pix, Δy + 0.25k max |err| {np.max(np.abs(vals[:, iy] + 0.25 * k)):.4f} pix, "
          f"dσ/(48(1+0.01k)) - 1 max {np.max(np.abs(vals[:, idn] / (48 * (1 + 0.01 * k)) - 1)):.4f}")
    # 2) export batch (csv) and a comparison (parquet)
    for f in glob.glob("/tmp/exp8*"): os.remove(f)
    paths = bw.export("/tmp/exp8_b"); print("2) batch export:", [os.path.basename(p) for p in paths])
    import pandas as pd; df = pd.read_csv(paths[0]); print("   csv columns:", list(df.columns)[:6], "... rows", len(df), "| dpos_x_pix equal:", np.allclose(df.dpos_x_pix, vals[:, ix]))
    app.cfg["export_format"] = "parquet"; app.paths = ["test_spot.fits"]
    def pick():
        d = [x for x in app.winfo_children() if isinstance(x, g.PairDialog)][0]; d._cb[0].current(0); d._cb[1].current(1); d._ok()
    app.after(300, pick); app.run(); wait(); rw = last(g.ResultWindow)
    paths = rw.export("/tmp/exp8_c"); print("   comparison export:", [os.path.basename(p) for p in paths], "| summary rows:", len(pd.read_parquet(paths[0])))
    # 3) busy state: a second job is ignored, buttons disabled
    app.paths = ["/tmp/drift.fits"]; app.cfg["batch_workers"] = 1; batch()
    print("3) during job: buttons", [b.instate(["disabled"]) for b in app._run_buttons], "| cancel enabled:", app.btn_cancel.instate(["!disabled"]))
    before = app._job; app.run(); print("   second run ignored:", app._job is before)
    # 4) cancel
    app.after(50, app.cancel_job); n_windows = len([w for w in app.winfo_children() if isinstance(w, g.BatchWindow)]); wait()
    print("4) after cancel:", app.status.get(), "| new batch window:", len([w for w in app.winfo_children() if isinstance(w, g.BatchWindow)]) > n_windows,
          "| buttons re-enabled:", [b.instate(["!disabled"]) for b in app._run_buttons])
    app.destroy()
