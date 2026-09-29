import os, sys, warnings; sys.path.insert(0, os.getcwd())
warnings.simplefilter("error")
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
app = g.App(); app.withdraw(); out = {}
print("GUI defaults:", app.v_nan.get(), app.v_taper.get())
app.paths = ["test_mef.fits"]
def pick(ia, ib):
    def f():
        dlg = [w for w in app.winfo_children() if isinstance(w, g.PairDialog)][0]
        dlg._cb[0].current(ia); dlg._cb[1].current(ib); dlg._ok()
    return f
for mode, taper in (("interp", 16.0), ("mask", 0.0)):
    app.v_nan.set(mode); app.v_taper.set(taper)
    app.after(300, pick(2, 4)); app.run(); wait()   # CUBE plane 0 vs MASKED
    rw = [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
    print("   Fourier title:", rw.figs[1][1]._suptitle.get_text())
    print(mode, "->", app.status.get())
# grab suptitle through the figure attached to the last canvas
import matplotlib.backends.backend_tkagg as bt
figs = [w.figure for w in bt.FigureCanvasTkAgg.__dict__.get("_instances", [])] if False else None
app.destroy(); print("ok")
