import os, sys, warnings; sys.path.insert(0, os.getcwd()); warnings.simplefilter("error")
import numpy as np
from matplotlib.backend_bases import MouseEvent
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_tables import SECTION
from lib.fitscmp_io import list_images, load_image
errors = []; g.messagebox.showerror = lambda *a, **k: errors.append(a[1])
app = g.App(); app.withdraw(); app.paths = ["test_spot.fits"]; app.cfg["saturation_level"] = 4095
last = lambda: [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
rows = lambda w: {r[0]: r[1] for r in w.table[1] if r[0] != SECTION}
def pick(ia, ib):
    def f():
        d = [x for x in app.winfo_children() if isinstance(x, g.PairDialog)][0]; d._cb[0].current(ia); d._cb[1].current(ib); d._ok()
    return f
# 1) without / with master dark: FARFIELD vs FARFIELD_HOT
app.after(300, pick(0, 4)); app.run(); wait(); r0 = rows(last())
app._set_dark("test_dark.fits"); print("1) dark loaded:", app.v_dark_lab.get(), "| subtract:", app.v_sub_dark.get())
app.after(300, pick(0, 4)); app.run(); wait(); w = last(); r1 = rows(w)
print(f"   B d_sigma without dark {r0['Width dσ [pix]'][1]:.2f}, with dark {r1['Width dσ [pix]'][1]:.2f} (truth 80); B background {r1['Background'][1]:.3f}")
print("   header:", w.header.splitlines()[-1][-60:])
# 2) saturation evaluated on raw values with dark subtraction
app.v_single.set(True)
def sel(i):
    def f():
        d = [x for x in app.winfo_children() if isinstance(x, g.SelectDialog)][0]; d._lb.selection_clear(0, "end"); d._lb.selection_set(i); d._ok()
    return f
app.after(300, sel(3)); app.run(); wait(); raw_sat = int((load_image(list_images("test_spot.fits")[3]) >= 4095).sum())
print("2) SATURATED with dark: saturated pixels", rows(last())["Saturated pixels"][0], "= raw count", raw_sat)
# 3) hot-pixel filter without dark
app.clear_dark(); app.v_hot.set(True); app.after(300, sel(4)); app.run(); wait()
print("3) hot filter: corrected", rows(last())["Hot pixels corrected"], "| d_sigma", round(rows(last())["Width dσ [pix]"][0], 2))
# 4) dark with wrong shape -> error, no window
app.v_hot.set(False); app._set_dark("test_dark.fits"); app.paths = ["test_mef.fits"]; n = len(errors)
app.after(300, sel(0)); app.run(); wait(); print("4) shape mismatch error:", errors[n:][:1])
# 5) readout on the comparison window
app.clear_dark(); app.paths = ["test_spot.fits"]; app.v_single.set(False); app.v_dy.set(1.0)
app.after(300, pick(0, 1)); app.run(); wait(); w = last()
fig = w.figs[0][1]; axA = [a for a in fig.axes if a.images][0]; canvas = fig.canvas
def event(ax, x, y, name="motion_notify_event", button=None):
    X, Y = ax.transData.transform((x, y)); return MouseEvent(name, ax.figure.canvas, X, Y, button=button)
txt = w.readout_text(event(axA, 300.2, 239.8)); print("5) readout at (300, 240):", txt)
a_raw = load_image(list_images("test_spot.fits")[0]); print("   A value from file:", a_raw[240, 300])
w._on_click(event(axA, 300, 240, "button_press_event", 1)); print("   pins after click:", len(w.pins))
canvas.toolbar.zoom(); w._on_click(event(axA, 310, 240, "button_press_event", 1)); print("   pins after click in zoom mode:", len(w.pins)); canvas.toolbar.zoom()
fax = [a for a in w.figs[1][1].axes if a.images][0]; print("   Fourier panel readout:", w.readout_text(event(fax, 0.1, 0.1)))
app.destroy()
