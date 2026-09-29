import os, sys, warnings; sys.path.insert(0, os.getcwd()); warnings.simplefilter("error")
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
from lib.fitscmp_tables import SECTION
app = g.App(); app.withdraw(); app.paths = ["test_spot.fits"]
last = lambda: [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
def pick(ia, ib):
    def f():
        d = [x for x in app.winfo_children() if isinstance(x, g.PairDialog)][0]; d._cb[0].current(ia); d._cb[1].current(ib); d._ok()
    return f
def rows(w): return {r[0]: r[1] for r in w.table[1] if r[0] != SECTION}
# ISO, auto registration on: deltas must be those of the raw images
app.v_auto.set(True); app.after(300, pick(0, 1)); app.run(); wait(); w = last(); r = rows(w)
print("columns:", w.table[0])
for k in ("Centroid x [pix]", "Centroid y [pix]", "Width dσx [pix]", "Width dσy [pix]", "Azimuth φ [deg]", "Ellipticity", "Total flux", "Round (ellipticity > 0.87)"):
    print(f"  {k:28s}", [_fmt(v) for v in r[k]])
print("  auto shift fields (dy, dx):", app.v_dy.get(), app.v_dx.get(), "| expected ≈ (-ΔY, -ΔX)")
print("  saturation row present without level:", "Saturated pixels" in r)
# Rainer method + saturation level from config
app.v_auto.set(False); app.reset_shift(); app.v_spot.set("Rainer (p10)"); app.cfg["saturation_level"] = 4095
app.v_single.set(True)
def sel():
    d = [x for x in app.winfo_children() if isinstance(x, g.SelectDialog)][0]; d._lb.selection_clear(0, "end"); d._lb.selection_set(0); d._lb.selection_set(3); d._ok()
app.after(300, sel); app.run(); wait(); w = last(); r = rows(w)
print("Rainer single view columns:", w.table[0], "| rows include FWHM:", "FWHM x [pix]" in r, "| d_sigma rows:", "Width dσ [pix]" in r)
print("  Saturated pixels:", r["Saturated pixels"], "| Centroid x:", [_fmt(v) for v in r["Centroid x [pix]"]])
print("  note starts:", w.table[2][:70])
app.destroy()
