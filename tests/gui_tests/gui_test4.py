import os, sys, warnings; sys.path.insert(0, os.getcwd()); warnings.simplefilter("error")
import numpy as np
from astropy.io import fits
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
from lib.fitscmp_analysis import image_stats
warned = []; g.messagebox.showwarning = lambda *a, **k: warned.append(a[1])
app = g.App(); app.withdraw()
last = lambda: [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
def images(win): return [ax.images[0].get_array() for f in (fig for _, fig in win.figs) for ax in f.axes if ax.images]
def select(indices):
    def f():
        dlg = [w for w in app.winfo_children() if isinstance(w, g.SelectDialog)][0]
        dlg._lb.selection_clear(0, "end")
        for i in indices: dlg._lb.selection_set(i)
        dlg._ok()
        if dlg.winfo_exists(): dlg.destroy()   # after a rejected selection
    return f
raw_single = fits.getdata("test_single.fits").astype(float)
# 1) a file with one image -> automatic single view, raw values
app.paths = ["test_single.fits"]; app.run(); wait(); w = last()
im = images(w); print("1) auto single:", w.title(), "| panels:", len(im), "| raw data shown:", np.array_equal(im[0], raw_single),
      "| status:", app.status.get())
tab = dict((lab, v) for lab, v in w.table[1]); print("   table max/sum vs numpy:", tab["Max"][0] == raw_single.max(), np.isclose(tab["Sum"][0], raw_single.sum()))
# 2) single mode with a multi-image file, button label
app.v_single.set(True); print("2) button label:", app.btn_run.cget("text"))
app.paths = ["test_mef.fits"]; app.after(300, select([0, 3, 4])); app.run(); wait(); w = last()
im = images(w); mef = [fits.getdata("test_mef.fits", 1), fits.getdata("test_mef.fits", 3)[1], fits.getdata("test_mef.fits", 4)]
print("   panels:", len(im), "| raw data:", all(np.array_equal(np.asarray(x), np.asarray(y, float), equal_nan=True) for x, y in zip(im, mef)),
      "| NaN count MASKED:", dict(w.table[1])["NaN pixels"][2], "=", int(np.isnan(mef[2]).sum()))
# 3) too many images -> warning, nothing shown
n_before = len([x for x in app.winfo_children() if isinstance(x, g.ResultWindow)])
app.after(300, select([0, 1, 2, 3, 4])); app.run(); wait()
print("3) 5 selected -> warning:", warned[-1:], "| new window:", len([x for x in app.winfo_children() if isinstance(x, g.ResultWindow)]) > n_before)
# 4) comparison mode keeps raw values in A and B panels with normalization
app.v_single.set(False); print("4) button label:", app.btn_run.cget("text")); app.v_norm.set("mean")
def pick():
    dlg = [x for x in app.winfo_children() if isinstance(x, g.PairDialog)][0]; dlg._cb[0].current(2); dlg._cb[1].current(0); dlg._ok()
app.after(300, pick); app.run(); wait(); w = last(); im = images(w)
a_raw = fits.getdata("test_mef.fits", 3)[0].astype(float); b_raw = fits.getdata("test_mef.fits", 1).astype(float)
print("   A raw:", np.array_equal(im[0], a_raw), "| B raw:", np.array_equal(im[1], b_raw), "| A-B normalized:", np.allclose(im[2], a_raw / a_raw.mean() - b_raw / b_raw.mean()),
      "| columns:", w.table[0], "| titles:", [ax.get_title() for ax in w.figs[0][1].axes if ax.images][:3])
app.destroy()
