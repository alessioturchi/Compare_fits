import os, sys, warnings; sys.path.insert(0, os.getcwd()); warnings.simplefilter("error")
import numpy as np
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
from lib.fitscmp_tables import SECTION

app = g.App(); app.withdraw(); app.paths = ["test_spot.fits"]; app.cfg["saturation_level"] = 4095; app.cfg["pixel_scale"] = 4.8  # as if set in fitscmp.yaml
last = lambda: [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
rows = lambda w: {r[0]: r[1] for r in w.table[1] if r[0] != SECTION}
def pick(ia, ib):
    def f():
        d = [x for x in app.winfo_children() if isinstance(x, g.PairDialog)][0]; d._cb[0].current(ia); d._cb[1].current(ib); d._ok()
    return f
def overlay_info(w):
    axes = [ax for ax in w.figs[0][1].axes if ax.images][:2]
    return [(len(ax.patches), [tuple(np.round(l.get_xydata()[0], 2)) for l in ax.lines if l.get_label() in ("centroid", "fit center")]) for ax in axes]
# 1) Gaussian fit, µm off, auto registration: overlay of B follows the shifted image
app.v_spot.set("Gaussian fit"); app.v_auto.set(True); app.after(300, pick(0, 1)); app.run(); wait(); w = last(); r = rows(w)
print("1) fit rows:", [k for k in r if k.startswith(("Fit center", "Width dσx", "Fit residual", "Saturated"))])
print("   Fit center x [pix]:", [_fmt(v) for v in r["Fit center x [pix]"]], "| error:", [_fmt(v) for v in r["Fit center x error [pix]"]])
print("   shift applied (dy, dx):", app.v_dy.get(), app.v_dx.get(), "| overlay (patches, center) A, B:", overlay_info(w))
# 2) µm units on ISO: positions and widths multiplied by 4.8, other rows unchanged
app.v_auto.set(False); app.reset_shift(); app.v_spot.set("ISO 11146")
app.after(300, pick(0, 1)); app.run(); wait(); r_pix = rows(last())
app.v_um.set(True); app.after(300, pick(0, 1)); app.run(); wait(); w = last(); r_um = rows(w)
ok = all(np.isclose(r_um[k.replace("[pix]", "[µm]")][i], 4.8 * r_pix[k][i]) for k in r_pix if "[pix]" in k for i in range(3) if isinstance(r_pix[k][i], float))
print("2) µm labels:", [k for k in r_um if "[µm]" in k][:4], "... | values = 4.8 x pix:", ok,
      "| flux ratio unchanged:", r_um["Total flux"][2] == r_pix["Total flux"][2], "| note:", "4.8 µm/pix" in w.table[2])
# 3) overlay off
app.v_overlay.set(False); app.after(300, pick(0, 1)); app.run(); wait(); print("3) overlay off -> patches/lines:", overlay_info(last()))
# 4) single view, Rainer, near field: FWHM ellipse, no area
app.v_overlay.set(True); app.v_single.set(True); app.v_spot.set("Rainer (p10)")
def sel():
    d = [x for x in app.winfo_children() if isinstance(x, g.SelectDialog)][0]; d._lb.selection_clear(0, "end"); d._lb.selection_set(2); d._ok()
app.after(300, sel); app.run(); wait(); w = last(); ax = [a for a in w.figs[0][1].axes if a.images][0]
print("4) Rainer overlay patches:", [type(p).__name__ for p in ax.patches], "| legend:", [t.get_text() for t in ax.get_legend().get_texts()])
app.destroy()
