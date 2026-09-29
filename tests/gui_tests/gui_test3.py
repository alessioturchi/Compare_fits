import os, sys, warnings; sys.path.insert(0, os.getcwd()); warnings.simplefilter("error")
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
app = g.App(); app.withdraw(); app.paths = ["test_mef.fits"]
def pick(ia, ib):
    def f():
        dlg = [w for w in app.winfo_children() if isinstance(w, g.PairDialog)][0]
        dlg._cb[0].current(ia); dlg._cb[1].current(ib); dlg._ok()
    return f
def run(ia, ib):
    app.after(300, pick(ia, ib)); app.run(); wait()
    rw = [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
    return rw.header.splitlines()[-1]
app.v_auto.set(True);  print("auto   CUBE0 vs SUBSHIFT:", run(2, 5), "| fields:", app.v_dy.get(), app.v_dx.get())
app.v_auto.set(False); app.v_dy.set(-1.25); print("manual edit dy       :", run(2, 5))
app.reset_shift();     print("after reset          :", run(2, 5))
app.v_auto.set(True);  print("auto   CUBE0 vs MASKED:", run(2, 4))
print("status:", app.status.get()[-60:])
app.destroy()
