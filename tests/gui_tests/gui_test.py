import os, sys, numpy as np
sys.path.insert(0, os.getcwd())
from astropy.io import fits
import fitscmp_gui as g
import time as _time
def wait():
    while app._job is not None:
        app.update(); _time.sleep(0.02)
from lib.fitscmp_windows import _fmt
fits.HDUList([fits.PrimaryHDU(np.random.rand(64, 80)), fits.ImageHDU(np.random.rand(64, 80), name="SECOND")]).writeto("/tmp/two.fits", overwrite=True)
app = g.App(); app.withdraw()
res = {}
def last_result():
    return [w for w in app.winfo_children() if isinstance(w, g.ResultWindow)][-1]
# 1) exactly two images, no swap / swap
app.paths = ["/tmp/two.fits"]
for sw in (False, True):
    app.v_swap.set(sw); app.run(); wait()
    lab = last_result().header
    res[f"two swap={sw}"] = lab.replace("\n", " | ")
# 2) dialog case: check default selection and accept
app.paths = ["test_mef.fits"]
def probe():
    dlg = [w for w in app.winfo_children() if isinstance(w, g.PairDialog)][0]
    res["dialog defaults (swap=True)"] = [cb.current() for cb in dlg._cb]
    dlg._ok()
app.v_swap.set(True); app.after(300, probe); app.run(); wait()
res["dialog result"] = app.status.get()
for k, v in res.items(): print(k, "->", v)
app.destroy()
