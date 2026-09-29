"""Tkinter GUI to view 2D images from FITS files and compare them (subtraction and Fourier analysis)."""
import os
import queue
import threading
import tkinter as tk
from tkinter import ttk, filedialog, messagebox

import yaml
import matplotlib
matplotlib.use("TkAgg")  # before importing the windows module (Tk canvases)

from lib.fitscmp_io import list_images
from lib.fitscmp_analysis import NORMALIZATIONS
from lib.fitscmp_calib import load_dark
from lib.fitscmp_export import EXPORT_FORMATS, export_batch, export_result
from lib.fitscmp_mask import NAN_MODES
from lib.fitscmp_pipeline import Cancelled, Settings, frc_crossing, run_batch, run_comparison, run_single
from lib.fitscmp_plots import plot_batch, plot_fourier, plot_single, plot_subtraction, spot_overlay
from lib.fitscmp_spot import SPOT_METHODS
from lib.fitscmp_tables import DEFINITIONS, batch_columns, batch_values, build_table
from lib.fitscmp_windows import BatchDialog, BatchWindow, HeaderWindow, PairDialog, ResultWindow, SelectDialog

CONFIG_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "fitscmp.yaml")
DEFAULTS = {"normalization": "none", "window": True, "radial_bins": 64,
            "clip_percentiles": [1.0, 99.0], "cmap": "viridis", "diff_cmap": "RdBu_r",
            "nan_mode": "interp", "mask_taper": 16.0, "interp_sigma": 1.0,
            "registration_upsample": 100, "spot_method": "iso", "saturation_level": None,
            "pixel_scale": None, "dark_frame": None, "hot_pixel_nsigma": 5.0,
            "export_format": "csv", "batch_workers": "auto",
            "expand_cubes": True, "start_dir": "."}
_METHOD_ID = {label: key for key, label in SPOT_METHODS.items()}  # GUI label -> method id


def load_config(path=CONFIG_PATH):
    """Load YAML config on top of defaults; invalid normalization/NaN mode fall back to defaults."""
    cfg = dict(DEFAULTS)
    if os.path.isfile(path):
        with open(path) as f:
            cfg.update(yaml.safe_load(f) or {})
    norm = str(cfg["normalization"]).lower()  # YAML 'off' would parse as False
    cfg["normalization"] = norm if norm in NORMALIZATIONS else "none"
    if str(cfg["nan_mode"]).lower() not in NAN_MODES:
        cfg["nan_mode"] = DEFAULTS["nan_mode"]
    if str(cfg["spot_method"]).lower() not in SPOT_METHODS:
        cfg["spot_method"] = DEFAULTS["spot_method"]
    if str(cfg["export_format"]).lower() not in EXPORT_FORMATS:
        cfg["export_format"] = DEFAULTS["export_format"]
    cfg["export_format"] = str(cfg["export_format"]).lower()
    return cfg


class App(tk.Tk):
    def __init__(self):
        super().__init__()
        self.title("FITS image comparison")
        self.cfg = load_config()
        self._last_dir, self.paths = self.cfg["start_dir"], []

        top = ttk.Frame(self, padding=5)
        top.pack(fill="both", expand=True)
        self.lb = tk.Listbox(top, selectmode="extended", width=90, height=12)
        self.lb.pack(side="left", fill="both", expand=True)
        btns = ttk.Frame(top)
        btns.pack(side="left", fill="y", padx=5)
        buttons = {}
        for text, cmd in (("Add files...", self.add_files), ("Remove", self.remove),
                          ("Clear", self.clear), ("Inspect", self.inspect),
                          ("Compare", self.run), ("Batch...", self.batch)):
            buttons[text] = ttk.Button(btns, text=text, command=cmd)
            buttons[text].pack(fill="x", pady=2)
        self.btn_run = buttons["Compare"]  # label follows the view mode
        self._run_buttons, self._job = (self.btn_run, buttons["Batch..."]), None

        opt = ttk.LabelFrame(self, text="Options", padding=5)
        opt.pack(fill="x", padx=5, pady=5)
        self.v_norm = tk.StringVar(value=self.cfg["normalization"])
        self.v_win = tk.BooleanVar(value=bool(self.cfg["window"]))
        self.v_bins = tk.IntVar(value=int(self.cfg["radial_bins"]))
        self.v_swap = tk.BooleanVar(value=False)
        self.v_single = tk.BooleanVar(value=False)
        self.v_single.trace_add("write", lambda *_: self.btn_run.configure(
            text="Show" if self.v_single.get() else "Compare"))
        ttk.Checkbutton(opt, text="Single images", variable=self.v_single).pack(side="left", padx=(0, 10))
        ttk.Label(opt, text="Normalization").pack(side="left")
        ttk.Combobox(opt, textvariable=self.v_norm, values=NORMALIZATIONS,
                     state="readonly", width=8).pack(side="left", padx=5)
        ttk.Checkbutton(opt, text="Hann window", variable=self.v_win).pack(side="left", padx=10)
        ttk.Label(opt, text="Radial bins").pack(side="left")
        ttk.Spinbox(opt, from_=8, to=1024, textvariable=self.v_bins, width=6).pack(side="left", padx=5)
        ttk.Checkbutton(opt, text="Swap A/B", variable=self.v_swap).pack(side="left", padx=10)

        cal = ttk.LabelFrame(self, text="Calibration (applied to all images)", padding=5)
        cal.pack(fill="x", padx=5, pady=(0, 5))
        self.dark, self.v_dark_lab = None, tk.StringVar(value="no master dark")
        self.v_sub_dark, self.v_hot = tk.BooleanVar(value=False), tk.BooleanVar(value=False)
        self._chk_dark = ttk.Checkbutton(cal, text="Subtract dark", variable=self.v_sub_dark)
        self._chk_dark.pack(side="left")
        self._chk_dark.state(["disabled"])
        ttk.Button(cal, text="Master dark...", command=self.select_dark).pack(side="left", padx=5)
        ttk.Button(cal, text="Clear", command=self.clear_dark).pack(side="left")
        ttk.Label(cal, textvariable=self.v_dark_lab).pack(side="left", padx=5)
        ttk.Checkbutton(cal, text=f"Hot-pixel filter ({self.cfg['hot_pixel_nsigma']:g} \u03c3)",
                        variable=self.v_hot).pack(side="left", padx=(15, 0))
        if self.cfg["dark_frame"]:
            self._set_dark(self.cfg["dark_frame"])

        msk = ttk.LabelFrame(self, text="Masked pixels (FFT)", padding=5)
        msk.pack(fill="x", padx=5, pady=(0, 5))
        self.v_nan = tk.StringVar(value=str(self.cfg["nan_mode"]).lower())
        self.v_taper = tk.DoubleVar(value=float(self.cfg["mask_taper"]))
        ttk.Label(msk, text="Mode").pack(side="left")
        ttk.Combobox(msk, textvariable=self.v_nan, values=NAN_MODES,
                     state="readonly", width=8).pack(side="left", padx=5)
        ttk.Label(msk, text="Taper [pix]").pack(side="left", padx=(10, 0))
        ttk.Spinbox(msk, from_=0, to=64, increment=1, textvariable=self.v_taper,
                    width=6).pack(side="left", padx=5)

        spot = ttk.LabelFrame(self, text="Spot analysis", padding=5)
        spot.pack(fill="x", padx=5, pady=(0, 5))
        self.v_spot = tk.StringVar(value=SPOT_METHODS[str(self.cfg["spot_method"]).lower()])
        ttk.Label(spot, text="Method").pack(side="left")
        ttk.Combobox(spot, textvariable=self.v_spot, values=list(SPOT_METHODS.values()),
                     state="readonly", width=14).pack(side="left", padx=5)
        self.v_overlay, self.v_um = tk.BooleanVar(value=True), tk.BooleanVar(value=False)
        ttk.Checkbutton(spot, text="Overlay", variable=self.v_overlay).pack(side="left", padx=(10, 0))
        scale = self.cfg["pixel_scale"]
        um = ttk.Checkbutton(spot, text="Units \u00b5m" + (f" ({scale:g} \u00b5m/pix)" if scale else " (pixel_scale not set)"),
                             variable=self.v_um)
        um.pack(side="left", padx=10)
        if not scale:
            um.state(["disabled"])
        sat = self.cfg["saturation_level"]
        ttk.Label(spot, text=f"Saturation level: {sat if sat is not None else 'not set'} (fitscmp.yaml)").pack(
            side="left", padx=10)

        reg = ttk.LabelFrame(self, text="Registration (shift applied to B)", padding=5)
        reg.pack(fill="x", padx=5, pady=(0, 5))
        self.v_dy, self.v_dx = tk.DoubleVar(value=0.0), tk.DoubleVar(value=0.0)
        self.v_auto = tk.BooleanVar(value=False)
        for text, var in (("dy [pix]", self.v_dy), ("dx [pix]", self.v_dx)):
            ttk.Label(reg, text=text).pack(side="left")
            ttk.Entry(reg, textvariable=var, width=8).pack(side="left", padx=(5, 10))
        ttk.Checkbutton(reg, text="Auto (estimate per comparison)", variable=self.v_auto).pack(side="left")
        ttk.Button(reg, text="Reset", command=self.reset_shift).pack(side="left", padx=10)

        bottom = ttk.Frame(self)
        bottom.pack(fill="x", padx=5, pady=(0, 5))
        self.status = tk.StringVar(value="Add FITS files. Inspect/Compare/Show/Batch act on the selected files (all if none).")
        self.btn_cancel = ttk.Button(bottom, text="Cancel", command=self.cancel_job)
        self.btn_cancel.pack(side="right")
        self.btn_cancel.state(["disabled"])
        self.progress = ttk.Progressbar(bottom, length=180)
        self.progress.pack(side="right", padx=5)
        ttk.Label(bottom, textvariable=self.status, anchor="w").pack(side="left", fill="x", expand=True)

    # ---- file list
    def add_files(self):
        new = filedialog.askopenfilenames(
            parent=self, initialdir=self._last_dir, title="Select FITS files",
            filetypes=[("FITS", "*.fits *.fit *.fts *.fz *.gz"), ("All files", "*")])
        for p in new:
            if p not in self.paths:
                self.paths.append(p)
                self.lb.insert("end", p)
        if new:
            self._last_dir = os.path.dirname(new[0])
        self.status.set(f"{len(self.paths)} file(s) in list.")

    def remove(self):
        for i in reversed(self.lb.curselection()):
            self.lb.delete(i)
            self.paths.pop(i)

    def clear(self):
        self.lb.delete(0, "end")
        self.paths.clear()

    def _targets(self):
        sel = self.lb.curselection()
        return [self.paths[i] for i in sel] if sel else list(self.paths)

    # ---- actions
    # ---- calibration
    def _set_dark(self, path):
        try:
            self.dark, label = load_dark(path)
        except Exception as e:
            messagebox.showerror("Master dark", f"{path}\n{e}", parent=self)
            return
        self.v_dark_lab.set(label)
        self._chk_dark.state(["!disabled"])
        self.v_sub_dark.set(True)

    def select_dark(self):
        path = filedialog.askopenfilename(parent=self, initialdir=self._last_dir, title="Select master dark",
                                          filetypes=[("FITS", "*.fits *.fit *.fts *.fz *.gz"), ("All files", "*")])
        if path:
            self._set_dark(path)

    def clear_dark(self):
        self.dark = None
        self.v_dark_lab.set("no master dark")
        self.v_sub_dark.set(False)
        self._chk_dark.state(["disabled"])

    # ---- settings and background jobs
    def settings(self):
        """Snapshot of the processing options, read in the GUI thread before a job starts."""
        dark = self.dark if self.v_sub_dark.get() else None
        return Settings(normalization=self.v_norm.get(), window=self.v_win.get(), radial_bins=int(self.v_bins.get()),
                        nan_mode=self.v_nan.get(), taper=float(self.v_taper.get()),
                        interp_sigma=float(self.cfg["interp_sigma"]), spot_method=_METHOD_ID[self.v_spot.get()],
                        saturation=self.cfg["saturation_level"], dark=dark,
                        dark_label=self.v_dark_lab.get() if dark is not None else "",
                        hot_filter=self.v_hot.get(), hot_nsigma=float(self.cfg["hot_pixel_nsigma"]),
                        auto_shift=self.v_auto.get(), dy=float(self.v_dy.get()), dx=float(self.v_dx.get()),
                        upsample=int(self.cfg["registration_upsample"]))

    def _scale(self):
        return self.cfg["pixel_scale"] if self.v_um.get() else None

    def _workers(self):
        """Number of batch processes: 'auto' = min(4, CPUs - 1), at least 1."""
        w = self.cfg["batch_workers"]
        return max(1, min(4, (os.cpu_count() or 1) - 1)) if str(w).lower() == "auto" else max(1, int(w))

    def _start(self, label, job, on_done, determinate=False):
        """Run job(progress, cancelled) in a worker thread; on_done(result) runs in the GUI thread.

        Tk is only used from the GUI thread: the worker communicates through a queue polled with after().
        """
        if self._job is not None:  # one job at a time
            return
        self._job = {"queue": queue.Queue(), "cancel": threading.Event()}
        q, cancel = self._job["queue"], self._job["cancel"]
        for b in self._run_buttons:
            b.state(["disabled"])
        self.btn_cancel.state(["!disabled"])
        self.progress.configure(mode="determinate" if determinate else "indeterminate", value=0, maximum=1)
        if not determinate:
            self.progress.start(15)
        self.status.set(label + "...")

        def work():
            try:
                q.put(("done", job(lambda done, total: q.put(("progress", done, total)), cancel.is_set)))
            except Cancelled:
                q.put(("cancelled",))
            except Exception as e:  # reported in the GUI thread
                q.put(("error", e))
        threading.Thread(target=work, daemon=True).start()
        self.after(100, self._poll, label, on_done)

    def _poll(self, label, on_done):
        q, cancel = self._job["queue"], self._job["cancel"]
        while True:
            try:
                msg = q.get_nowait()
            except queue.Empty:
                self.after(100, self._poll, label, on_done)
                return
            if msg[0] != "progress":
                break
            self.progress.configure(value=msg[1], maximum=msg[2])
            self.status.set(f"{label}: {msg[1]}/{msg[2]}")
        self._finish()
        if msg[0] == "error":
            self.status.set(f"{label}: error")
            messagebox.showerror(label, str(msg[1]), parent=self)
        elif msg[0] == "cancelled" or cancel.is_set():  # a cancelled single job is discarded
            self.status.set(f"{label}: cancelled")
        else:
            try:
                on_done(msg[1])
            except Exception as e:
                messagebox.showerror(label, str(e), parent=self)

    def _finish(self):
        self.progress.stop()
        self.progress.configure(value=0)
        for b in self._run_buttons:
            b.state(["!disabled"])
        self.btn_cancel.state(["disabled"])
        self._job = None

    def cancel_job(self):
        if self._job is not None:
            self._job["cancel"].set()
            self.status.set("Cancelling...")

    # ---- actions
    def reset_shift(self):
        self.v_dy.set(0.0)
        self.v_dx.set(0.0)

    def inspect(self):
        paths = self._targets()
        if not paths:
            messagebox.showinfo("Inspect", "No FITS file in the list.", parent=self)
            return
        HeaderWindow(self, paths)

    def _collect(self):
        """All 2D images of the target files, or None after an error message."""
        refs = []
        for p in self._targets():
            try:
                refs += list_images(p, self.cfg["expand_cubes"])
            except Exception as e:
                messagebox.showerror("Read error", f"{p}\n{e}", parent=self)
                return None
        if not refs:
            messagebox.showerror("No image", "No 2D image found in the selected files.", parent=self)
            return None
        return refs

    def run(self):
        """Single-image view if selected or if only one image is available, comparison otherwise."""
        if self._job is not None:  # busy: no dialog while a job runs
            return
        refs = self._collect()
        if refs is None:
            return
        if self.v_single.get() or len(refs) == 1:
            self.show_single(refs)
        else:
            self.compare(refs)

    def _table(self, imgs, names, res, s, delta=False):
        """Side-panel table (columns, rows, note) of a single-view or comparison result."""
        extra = [("Hot pixels corrected", res["n_hot"])] if s.hot_filter else []
        rows, note = build_table(imgs, names, res["spots"], s.spot_method, s.saturation, delta, self._scale(), extra)
        return (names + (["B \u2212 A"] if delta else []), rows, note)

    def show_single(self, refs):
        sel = refs if len(refs) == 1 else SelectDialog(self, refs).result
        if not sel:
            return
        s = self.settings()
        auto = " (only one image available)" if len(refs) == 1 and not self.v_single.get() else ""
        self._start("Loading images", lambda progress, cancelled: run_single(sel, s),
                    lambda res: self._show_single_result(sel, s, res, auto))

    def _show_single_result(self, sel, s, res, auto):
        imgs, names = res["imgs"], [f"#{k}" for k in range(len(sel))]
        table = self._table(imgs, names, res, s)
        spots = [m for m, _ in res["spots"]]
        overlays = [spot_overlay(m) for m in spots] if self.v_overlay.get() else None
        fig = plot_single(imgs, [f"#{k}  {r.label()}" for k, r in enumerate(sel)], self.cfg, overlays)
        header = ("\n".join(f"#{k}: {r.label()}" for k, r in enumerate(sel))
                  + f"\nValues in original units; calibration: {s.calibration_text()}")
        meta = {"mode": "single", "images": {n: r.label() for n, r in zip(names, sel)}, "settings": s.describe(),
                "pixel_scale_um": self._scale(), "hot_pixels_corrected": res["n_hot"],
                "spots": {n: m for n, m in zip(names, spots)}}
        fmt = self.cfg["export_format"]
        ResultWindow(self, "FITS images", header, [("Images", fig)], table,
                     readout=list(zip(names, imgs)), scale=self._scale(),
                     export=lambda base: export_result(base, table, meta, None, fmt))
        self.status.set(f"Shown {len(sel)} image(s){auto}")

    def compare(self, refs):
        swap = self.v_swap.get()
        if len(refs) == 2:
            pair = tuple(refs[::-1]) if swap else tuple(refs)
        else:
            pair = PairDialog(self, refs, swap).result  # swap only sets the default selection
        if not pair:
            return
        s = self.settings()
        self._start("Comparing", lambda progress, cancelled: run_comparison(pair, s),
                    lambda res: self._show_comparison(pair, s, res))

    def _show_comparison(self, pair, s, res):
        dy, dx, fc = res["dy"], res["dx"], res["fc"]
        if s.auto_shift:  # show the estimate (3 decimals, the applied value) in the fields
            self.v_dy.set(dy)
            self.v_dx.set(dx)
        table = self._table([res["a"], res["b"]], ["A", "B"], res, s, delta=True)  # B before the shift
        spots = [m for m, _ in res["spots"]]
        # overlays follow the displayed images: B is shown after the shift
        overlays = (spot_overlay(spots[0]), spot_overlay(spots[1], dx, dy)) if self.v_overlay.get() else (None, None)
        info = (f"B shifted by (dy, dx) = ({dy:+.3f}, {dx:+.3f}) pix" + (" [auto]" if s.auto_shift else "")
                + f"; calibration: {s.calibration_text()}")
        note = (f"NaN mode: {s.nan_mode}, taper {s.taper:g} pix, "
                f"pixels excluded from FFT: {100 * fc['excluded_frac']:.2f}%")
        figs = [("Subtraction", plot_subtraction(res["a"], res["b_shifted"], res["d"], res["stats"], self.cfg,
                                                 s.normalization, overlays)),
                ("Fourier", plot_fourier(fc, self.cfg, note))]
        header = f"A: {pair[0].label()}\nB: {pair[1].label()}\n{info}"
        readout = [("A", res["a"]), ("B" if dy == dx == 0 else "B(shifted)", res["b_shifted"]),
                   ("A-B" if s.normalization == "none" else f"A-B({s.normalization})", res["d"])]
        meta = {"mode": "comparison", "images": {"A": pair[0].label(), "B": pair[1].label()}, "settings": s.describe(),
                "applied_shift_pix": {"dy": dy, "dx": dx}, "pixel_scale_um": self._scale(),
                "hot_pixels_corrected": res["n_hot"], "difference_stats": res["stats"],
                "fourier": {"excluded_frac": fc["excluded_frac"], "frc_f50": frc_crossing(fc["f"], fc["frc"])},
                "spots": {"A": spots[0], "B": spots[1]}}
        fmt = self.cfg["export_format"]
        ResultWindow(self, "FITS comparison", header, figs, table, readout=readout, scale=self._scale(),
                     export=lambda base: export_result(base, table, meta, fc, fmt))
        self.status.set(f"Compared {pair[0].label()}  vs  {pair[1].label()}; {info}")

    def batch(self):
        """Compare all images of the target files with a reference chosen in a dialog."""
        if self._job is not None:  # busy: no dialog while a job runs
            return
        refs = self._collect()
        if refs is None:
            return
        if len(refs) < 2:
            messagebox.showerror("Batch", "At least two 2D images are needed.", parent=self)
            return
        k = BatchDialog(self, refs).result
        if k is None:
            return
        ref, others, s, workers = refs[k], [r for i, r in enumerate(refs) if i != k], self.settings(), self._workers()
        self._start(f"Batch of {len(others)} image(s)",
                    lambda progress, cancelled: run_batch(ref, others, s, workers, progress, cancelled),
                    lambda rows: self._show_batch(ref, s, rows, workers), determinate=True)

    def _show_batch(self, ref, s, rows, workers):
        scale = self._scale()
        cols = batch_columns(s.spot_method, s.saturation, s.hot_filter, scale)
        values = batch_values(rows, cols, scale)
        reg = "automatic per image" if s.auto_shift else f"fixed ({s.dy:+.3f}, {s.dx:+.3f}) pix"
        header = (f"Reference (A, index 0): {ref.label()}\n{len(rows) - 1} image(s) compared; "
                  f"spot: {SPOT_METHODS[s.spot_method]}; calibration: {s.calibration_text()}; registration: {reg}")
        note = (f"\u0394x, \u0394y and flux ratio: spot of each image minus (over) the reference, as recorded "
                f"(before registration). Shift: applied registration. rms_rel, Pearson: difference after "
                f"registration and normalization ({s.normalization}). FRC 0.5: first spatial frequency where the "
                f"Fourier Ring Correlation with the reference drops below 0.5. " + DEFINITIONS[s.spot_method])
        meta = {"mode": "batch", "reference": ref.label(), "settings": s.describe(), "pixel_scale_um": scale,
                "workers": workers, "rows": rows}
        fmt = self.cfg["export_format"]
        BatchWindow(self, header, [c[1] for c in cols], values, plot_batch(rows, s.spot_method, scale), note,
                    export=lambda base: export_batch(base, [c[2] for c in cols], values, meta, fmt))
        self.status.set(f"Batch done: {len(rows) - 1} image(s) compared with {ref.label()}")


if __name__ == "__main__":
    App().mainloop()
