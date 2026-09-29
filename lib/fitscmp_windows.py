"""Tkinter windows of fitscmp: selection dialogs, header browser, result and batch windows."""
import os
import tkinter as tk
from tkinter import filedialog, messagebox, ttk
from types import SimpleNamespace

import numpy as np
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk

from .fitscmp_io import read_headers
from .fitscmp_tables import SECTION

MAX_SINGLE = 4  # maximum number of images in the single-image view
MAX_PINS = 5    # pinned cursor readouts kept in the result window


def _refit(canvas):
    """Re-apply the widget size to the figure.

    Works around a Matplotlib TkAgg ordering issue on Linux screens whose DPI differs
    from 96: the device pixel ratio is set after the first resize, and the figure
    stays larger than a canvas that cannot grow (here, next to the side table).
    """
    tkw = canvas.get_tk_widget()
    if tkw.winfo_width() > 1:
        canvas.resize(SimpleNamespace(width=tkw.winfo_width(), height=tkw.winfo_height()))


def embed_figure(win, parent, fig):
    """Figure canvas with navigation toolbar packed in `parent`; returns the canvas."""
    canvas = FigureCanvasTkAgg(fig, master=parent)
    NavigationToolbar2Tk(canvas, parent).update()
    # small requested size: the canvas (and figure) adapt to the available space
    canvas.get_tk_widget().configure(width=600, height=400)
    canvas.get_tk_widget().pack(fill="both", expand=True)
    canvas.get_tk_widget().bind("<Map>", lambda _e: win.after_idle(_refit, canvas), add="+")
    canvas.draw()
    return canvas


def ask_export(win, export, initial="fitscmp"):
    """Ask a base file name and run export(base) -> list of written paths; report the result."""
    path = filedialog.asksaveasfilename(parent=win, title="Export results: base file name", initialfile=initial)
    if not path:
        return None
    base = os.path.splitext(path)[0] if os.path.splitext(path)[1].lower() in (".csv", ".parquet", ".json") else path
    try:
        paths = export(base)
    except Exception as e:
        messagebox.showerror("Export error", str(e), parent=win)
        return None
    messagebox.showinfo("Export", "Written:\n" + "\n".join(paths), parent=win)
    return paths


def _header_bar(win, header, export, initial):
    """Header text with an optional Export button on the right."""
    top = ttk.Frame(win)
    top.pack(fill="x", padx=5, pady=3)
    ttk.Label(top, text=header, justify="left").pack(side="left", anchor="w")
    if export:
        ttk.Button(top, text="Export...", command=lambda: ask_export(win, export, initial)).pack(side="right", anchor="n")


def _fmt(v):
    """Table cell text: numbers with 6 significant digits, booleans as yes/no, None and NaN as a dash."""
    if v is None or (isinstance(v, float) and v != v):
        return "\u2014"
    if isinstance(v, bool):
        return "yes" if v else "no"
    if isinstance(v, str):
        return v
    return f"{v:.6g}"



class PairDialog(tk.Toplevel):
    """Modal dialog to choose images A and B among more than two candidates."""

    def __init__(self, parent, refs, swap=False):
        super().__init__(parent)
        self.title("Select images to compare")
        self.transient(parent)
        self.result, self._refs, self._cb = None, refs, []
        labels = [f"#{i}  {r.label()}" for i, r in enumerate(refs)]
        width = min(max(len(s) for s in labels), 100)
        defaults = (1, 0) if swap else (0, 1)  # default A/B indices
        for row, (name, default) in enumerate(zip(("Image A", "Image B"), defaults)):
            ttk.Label(self, text=name).grid(row=row, column=0, padx=5, pady=5, sticky="w")
            cb = ttk.Combobox(self, values=labels, state="readonly", width=width)
            cb.current(default)
            cb.grid(row=row, column=1, padx=5, pady=5)
            self._cb.append(cb)
        fr = ttk.Frame(self)
        fr.grid(row=2, column=0, columnspan=2, pady=5)
        ttk.Button(fr, text="OK", command=self._ok).pack(side="left", padx=5)
        ttk.Button(fr, text="Cancel", command=self.destroy).pack(side="left", padx=5)
        self.grab_set()
        self.wait_window()

    def _ok(self):
        ia, ib = (cb.current() for cb in self._cb)
        if ia == ib:
            messagebox.showwarning("Selection", "A and B must be different images.", parent=self)
            return
        self.result = (self._refs[ia], self._refs[ib])
        self.destroy()


class SelectDialog(tk.Toplevel):
    """Modal dialog to choose 1 to MAX_SINGLE images for the single-image view."""

    def __init__(self, parent, refs):
        super().__init__(parent)
        self.title("Select images to show")
        self.transient(parent)
        self.result, self._refs = None, refs
        labels = [f"#{i}  {r.label()}" for i, r in enumerate(refs)]
        ttk.Label(self, text=f"Select 1 to {MAX_SINGLE} images (Ctrl/Shift for multiple)").pack(
            anchor="w", padx=5, pady=(5, 0))
        self._lb = tk.Listbox(self, selectmode="extended", exportselection=False,
                              width=min(max(len(s) for s in labels), 100), height=min(len(labels), 15))
        for s in labels:
            self._lb.insert("end", s)
        self._lb.selection_set(0)
        self._lb.pack(fill="both", expand=True, padx=5, pady=5)
        fr = ttk.Frame(self)
        fr.pack(pady=5)
        ttk.Button(fr, text="OK", command=self._ok).pack(side="left", padx=5)
        ttk.Button(fr, text="Cancel", command=self.destroy).pack(side="left", padx=5)
        self.grab_set()
        self.wait_window()

    def _ok(self):
        sel = self._lb.curselection()
        if not 1 <= len(sel) <= MAX_SINGLE:
            messagebox.showwarning("Selection", f"Select 1 to {MAX_SINGLE} images.", parent=self)
            return
        self.result = [self._refs[i] for i in sel]
        self.destroy()


class HeaderWindow(tk.Toplevel):
    """Header browser: one list entry per (file, HDU)."""

    def __init__(self, parent, paths):
        super().__init__(parent)
        self.title("FITS headers")
        self.geometry("1000x600")
        self._items = []
        for p in paths:
            try:
                self._items += [(f"{os.path.basename(p)} | {t}", h) for t, h in read_headers(p)]
            except Exception as e:
                self._items.append((f"{os.path.basename(p)} | ERROR", str(e)))
        pw = ttk.PanedWindow(self, orient="horizontal")
        pw.pack(fill="both", expand=True)
        self._lb = tk.Listbox(pw, width=45, exportselection=False)
        for name, _ in self._items:
            self._lb.insert("end", name)
        fr = ttk.Frame(pw)
        self._txt = tk.Text(fr, font=("Courier", 10), wrap="none")
        sy = ttk.Scrollbar(fr, orient="vertical", command=self._txt.yview)
        self._txt.configure(yscrollcommand=sy.set)
        sy.pack(side="right", fill="y")
        self._txt.pack(fill="both", expand=True)
        pw.add(self._lb, weight=1)
        pw.add(fr, weight=3)
        self._lb.bind("<<ListboxSelect>>", self._show)
        if self._items:
            self._lb.selection_set(0)
            self._show()

    def _show(self, _event=None):
        sel = self._lb.curselection()
        if not sel:
            return
        self._txt.configure(state="normal")
        self._txt.delete("1.0", "end")
        self._txt.insert("1.0", self._items[sel[0]][1])
        self._txt.configure(state="disabled")


class ResultWindow(tk.Toplevel):
    """Figures in tabs with matplotlib toolbars, an optional side table and a pixel readout.

    readout: [(name, array)] of the images as displayed; moving the cursor over an image
    panel shows the value of every array at that pixel, a left click pins the readout.
    scale: pixel scale [um/pix] to also show the position in um.
    """

    def __init__(self, parent, title, header, figs, table=None, readout=None, scale=None, export=None):
        super().__init__(parent)
        self.title(title)
        self.geometry("1500x900")
        self.figs, self.table, self.export, self.header = figs, table, export, header  # for inspection and tests
        self._readout, self._scale, self.pins = readout or [], scale, []
        _header_bar(self, header, export, "fitscmp")
        if self._readout:
            bar = ttk.Frame(self)
            bar.pack(fill="x", padx=5)
            self.v_read = tk.StringVar(value="Cursor on an image: pixel values in original units "
                                             "(left click to pin, when zoom/pan are off)")
            self.v_pins = tk.StringVar(value="")
            ttk.Label(bar, textvariable=self.v_read, font=("TkFixedFont", 9)).pack(anchor="w")
            ttk.Label(bar, textvariable=self.v_pins, font=("TkFixedFont", 9), foreground="#404040",
                      justify="left").pack(anchor="w")
            ttk.Button(bar, text="Clear pins", command=self.clear_pins).pack(anchor="w", pady=(0, 3))
        pw = ttk.PanedWindow(self, orient="horizontal")
        pw.pack(fill="both", expand=True)
        nb = ttk.Notebook(pw)
        pw.add(nb, weight=4)
        if table:
            pw.add(self._table(pw, *table), weight=1)
        for name, fig in figs:
            fr = ttk.Frame(nb)
            nb.add(fr, text=name)
            canvas = embed_figure(self, fr, fig)
            if self._readout:
                canvas.mpl_connect("motion_notify_event", self._on_move)
                canvas.mpl_connect("button_press_event", self._on_click)

    # ---- pixel readout
    def readout_text(self, event):
        """Readout string for a Matplotlib mouse event, or None outside image panels."""
        ax = event.inaxes
        if ax is None or not ax.images or event.xdata is None:
            return None
        l, r, b, t = ax.images[0].get_extent()
        shape = (int(round(t - b)), int(round(r - l)))
        if (l, b) != (-0.5, -0.5):  # pixel panels only (not the Fourier spectra)
            return None
        ix, iy = int(round(event.xdata)), int(round(event.ydata))
        if not (0 <= ix < shape[1] and 0 <= iy < shape[0]):
            return None
        pos = f"x={ix:4d} y={iy:4d}"
        if self._scale:
            pos += f" ({ix * self._scale:.1f}, {iy * self._scale:.1f} \u00b5m)"
        vals = [f"{name}={_fmt(float(arr[iy, ix])) if np.isfinite(arr[iy, ix]) else 'NaN'}"
                for name, arr in self._readout if arr.shape == shape]
        return "  ".join([pos] + vals)

    def _on_move(self, event):
        text = self.readout_text(event)
        if text:
            self.v_read.set(text)

    def _on_click(self, event):
        toolbar = event.canvas.toolbar
        if event.button != 1 or (toolbar is not None and toolbar.mode):
            return  # zoom/pan active: the click belongs to the toolbar
        text = self.readout_text(event)
        if text:
            self.pins = (self.pins + [text])[-MAX_PINS:]
            self.v_pins.set("\n".join(f"pin {k + 1}: {p}" for k, p in enumerate(self.pins)))

    def clear_pins(self):
        self.pins = []
        self.v_pins.set("")

    @staticmethod
    def _table(parent, columns, rows, note=""):
        """Read-only table (quantity + one column per image) with section rows and a note below."""
        fr = ttk.Frame(parent)
        ids = [f"c{k}" for k in range(len(columns))]  # "#N" ids are reserved by Treeview
        tv = ttk.Treeview(fr, columns=["q"] + ids, show="headings", height=len(rows))
        tv.heading("q", text="Quantity")
        tv.column("q", width=190, anchor="w", stretch=False)
        for cid, name in zip(ids, columns):
            tv.heading(cid, text=name)
            tv.column(cid, width=100, anchor="e")
        tv.tag_configure(SECTION, background="#e4e4e4", font=("TkDefaultFont", 9, "bold"))
        for row in rows:
            if row[0] == SECTION:
                tv.insert("", "end", values=[row[1]] + [""] * len(ids), tags=(SECTION,))
            else:
                tv.insert("", "end", values=[row[0]] + [_fmt(v) for v in row[1]])
        tv.pack(fill="both", expand=True, padx=5, pady=5)
        if note:
            msg = ttk.Label(fr, text=note, justify="left", wraplength=400, foreground="#404040")
            msg.pack(fill="x", padx=5, pady=(0, 5))
            fr.bind("<Configure>", lambda e: msg.configure(wraplength=max(e.width - 10, 100)), add="+")
        return fr


class BatchDialog(tk.Toplevel):
    """Modal dialog to choose the reference image of a batch comparison."""

    def __init__(self, parent, refs):
        super().__init__(parent)
        self.title("Batch comparison")
        self.transient(parent)
        self.result = None
        labels = [f"#{i}  {r.label()}" for i, r in enumerate(refs)]
        ttk.Label(self, text="Reference (A)").grid(row=0, column=0, padx=5, pady=5, sticky="w")
        self._cb = ttk.Combobox(self, values=labels, state="readonly", width=min(max(len(s) for s in labels), 100))
        self._cb.current(0)
        self._cb.grid(row=0, column=1, padx=5, pady=5)
        ttk.Label(self, text=f"The other {len(refs) - 1} image(s) are compared with the reference (B).").grid(
            row=1, column=0, columnspan=2, padx=5, sticky="w")
        fr = ttk.Frame(self)
        fr.grid(row=2, column=0, columnspan=2, pady=5)
        ttk.Button(fr, text="Run", command=self._ok).pack(side="left", padx=5)
        ttk.Button(fr, text="Cancel", command=self.destroy).pack(side="left", padx=5)
        self.grab_set()
        self.wait_window()

    def _ok(self):
        self.result = self._cb.current()
        self.destroy()


class BatchWindow(tk.Toplevel):
    """Batch results: table (one row per image) and trend figure, with export."""

    def __init__(self, parent, header, columns, values, fig, note="", export=None):
        super().__init__(parent)
        self.title("FITS batch comparison")
        self.geometry("1500x900")
        self.columns, self.values, self.fig, self.export, self.header = columns, values, fig, export, header  # for tests
        _header_bar(self, header, export, "fitscmp_batch")
        if note:  # packed before the notebook, at the bottom, so that it is never cut
            ttk.Label(self, text=note, justify="left", foreground="#404040", wraplength=1400).pack(
                side="bottom", fill="x", padx=5, pady=3)
        nb = ttk.Notebook(self)
        nb.pack(fill="both", expand=True)
        tab = ttk.Frame(nb)
        nb.add(tab, text="Table")
        ids = [f"c{k}" for k in range(len(columns))]  # "#N" ids are reserved by Treeview
        tv = ttk.Treeview(tab, columns=ids, show="headings")
        for cid, name in zip(ids, columns):
            tv.heading(cid, text=name)
            width = 320 if name == "Image" else max(100, 9 * len(name))  # header must stay readable
            tv.column(cid, width=width, anchor="w" if name in ("Image", "Error") else "e", stretch=name == "Error")
        for vals in values:
            tv.insert("", "end", values=[_fmt(v) if not isinstance(v, str) else v for v in vals])
        sx = ttk.Scrollbar(tab, orient="horizontal", command=tv.xview)
        sy = ttk.Scrollbar(tab, orient="vertical", command=tv.yview)
        tv.configure(xscrollcommand=sx.set, yscrollcommand=sy.set)
        sy.pack(side="right", fill="y")
        sx.pack(side="bottom", fill="x")
        tv.pack(fill="both", expand=True)
        fr = ttk.Frame(nb)
        nb.add(fr, text="Trends")
        embed_figure(self, fr, fig)
