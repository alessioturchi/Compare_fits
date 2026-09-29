"""Matplotlib figures for fitscmp (backend-independent: only matplotlib.figure is used)."""
import numpy as np
from matplotlib.figure import Figure
from matplotlib.patches import Ellipse, Polygon


def _clim(p, *imgs):
    """Common color limits from percentiles of the finite pixels of all images."""
    v = np.concatenate([i[np.isfinite(i)] for i in imgs])
    return tuple(np.percentile(v, p)) if v.size else (None, None)


def _sym(p_hi, img):
    """Symmetric color limits for signed maps."""
    s = _clim([p_hi], np.abs(img))[0]
    return (-s, s) if s else (None, None)


def spot_overlay(m, dx=0.0, dy=0.0):
    """Overlay geometry from spot metrics, optionally shifted by (dx, dy) pixels; None if no metrics."""
    if m is None:
        return None
    if m["method"] == "gauss":
        c, lab = (m["center_x"] + dx, m["center_y"] + dy), "fit center"
    else:
        c, lab = (m["centroid_x"] + dx, m["centroid_y"] + dy), "centroid"
    ov = {"center": c, "center_label": lab, "peak": (m["peak_x"] + dx, m["peak_y"] + dy)}
    if m["method"] == "rainer":  # FWHM ellipse along the frame axes
        ov["ellipse"], ov["ellipse_label"] = (m["fwhm_x"], m["fwhm_y"], 0.0), "FWHM"
    else:  # 2 sigma ellipse: full axes = d_sigma_x, d_sigma_y along the principal axes
        phi = 0.0 if m["phi_deg"] is None else m["phi_deg"]
        ov["ellipse"], ov["ellipse_label"] = (m["d_sigma_x"], m["d_sigma_y"], phi), "dσ ellipse"
    if m["method"] == "iso":  # integration area: rectangle of sides 3 x widths
        au, av, ap = m["area"]
        u, v = np.array([1, 1, -1, -1]) * 1.5 * au, np.array([1, -1, -1, 1]) * 1.5 * av
        ov["area"] = np.column_stack([c[0] + u * np.cos(ap) - v * np.sin(ap), c[1] + u * np.sin(ap) + v * np.cos(ap)])
        ov["area_label"] = "integration area"
    elif m["method"] == "gauss":  # fit box, pixel edges
        x0, x1, y0, y1 = m["fit_box"]
        ov["area"] = np.array([[x0, y0], [x1, y0], [x1, y1], [x0, y1]]) - 0.5 + [dx, dy]
        ov["area_label"] = "fit region"
    return ov


def _draw_overlay(ax, ov):
    """Center (+), peak position (x), ellipse and integration area on an image panel."""
    lim = ax.get_xlim(), ax.get_ylim()
    w, h, ang = ov["ellipse"]
    ax.add_patch(Ellipse(ov["center"], w, h, angle=ang, fill=False, ec="red", lw=1.2, label=ov["ellipse_label"]))
    if "area" in ov:
        ax.add_patch(Polygon(ov["area"], closed=True, fill=False, ec="white", ls="--", lw=1, label=ov["area_label"]))
    ax.plot(*ov["center"], "+", color="red", ms=12, mew=1.5, label=ov["center_label"])
    ax.plot(*ov["peak"], "x", color="white", ms=8, mew=1.5, label="peak position")
    ax.set_xlim(*lim[0])  # the overlay must not change the image limits
    ax.set_ylim(*lim[1])
    ax.legend(fontsize=7, loc="upper right", framealpha=0.6)


def _imshow(fig, ax, img, title, cmap, clim, extent=None, overlay=None):
    im = ax.imshow(img, origin="lower", cmap=cmap, vmin=clim[0], vmax=clim[1],
                   extent=extent, interpolation="nearest")
    ax.set_title(title, fontsize=9)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    if overlay:
        _draw_overlay(ax, overlay)


def plot_subtraction(a, b, d, stats, cfg, norm="none", overlays=(None, None)):
    """A and B in original units (with optional spot overlays); A - B and its statistics after `norm`."""
    fig = Figure(figsize=(11, 8), layout="constrained")
    ax1 = fig.add_subplot(2, 2, 1)
    ax2 = fig.add_subplot(2, 2, 2, sharex=ax1, sharey=ax1)
    ax3 = fig.add_subplot(2, 2, 3, sharex=ax1, sharey=ax1)
    ax4 = fig.add_subplot(2, 2, 4)
    p = cfg["clip_percentiles"]
    if norm == "none":  # same units: common color scale
        clim_a = clim_b = _clim(p, a, b)
    else:  # different flux levels: independent scales
        clim_a, clim_b = _clim(p, a), _clim(p, b)
    _imshow(fig, ax1, a, "A (original units)", cfg["cmap"], clim_a, overlay=overlays[0])
    _imshow(fig, ax2, b, "B (original units)", cfg["cmap"], clim_b, overlay=overlays[1])
    suffix = "" if norm == "none" else f" (normalized: {norm})"
    _imshow(fig, ax3, d, "A - B" + suffix, cfg["diff_cmap"], _sym(p[1], d))
    ax4.hist(d[np.isfinite(d)], bins=200, log=True, color="gray")
    ax4.set_title("Histogram of A - B" + suffix, fontsize=9)
    txt = "\n".join(f"{k:8s} {v:.4g}" for k, v in stats.items())
    ax4.text(0.02, 0.98, txt, transform=ax4.transAxes, va="top", fontsize=8,
             family="monospace", bbox=dict(fc="white", alpha=0.8))
    return fig


def plot_fourier(fc, cfg, note=""):
    fig = Figure(figsize=(14, 8), layout="constrained")
    if note:
        fig.suptitle(note, fontsize=9)
    ax = [fig.add_subplot(2, 3, 1)]
    ax += [fig.add_subplot(2, 3, k, sharex=ax[0], sharey=ax[0]) for k in (2, 3)]
    ax += [fig.add_subplot(2, 3, k) for k in (4, 5, 6)]
    p, fx, fy = cfg["clip_percentiles"], fc["fx"], fc["fy"]
    hx, hy = 0.5 / len(fx), 0.5 / len(fy)  # half frequency step: extent at pixel edges
    ext = [fx[0] - hx, fx[-1] + hx, fy[0] - hy, fy[-1] + hy]
    clim = _clim(p, fc["logps_a"], fc["logps_b"])
    _imshow(fig, ax[0], fc["logps_a"], "log10 P(A)", cfg["cmap"], clim, ext)
    _imshow(fig, ax[1], fc["logps_b"], "log10 P(B)", cfg["cmap"], clim, ext)
    _imshow(fig, ax[2], fc["log_ratio"], "log10 P(A)/P(B)", cfg["diff_cmap"],
            _sym(p[1], fc["log_ratio"]), ext)
    for a in ax[:3]:
        a.set_xlabel("fx [cycles/pix]")
    ax[0].set_ylabel("fy [cycles/pix]")

    f = fc["f"]
    for key, lab in (("prof_a", "A"), ("prof_b", "B"), ("prof_d", "A - B")):
        ax[3].loglog(f, fc[key], label=lab)
    ax[3].set(title="Radial power spectrum", xlabel="f [cycles/pix]", ylabel="power")
    ax[3].legend(fontsize=8)
    ax[4].loglog(f, fc["prof_a"] / fc["prof_b"])
    ax[4].axhline(1, ls=":", c="k")
    ax[4].set(title="Radial ratio P(A)/P(B)", xlabel="f [cycles/pix]")
    ax[5].semilogx(f, fc["frc"])
    ax[5].axhline(0, ls=":", c="k")
    ax[5].set(title="Fourier Ring Correlation", xlabel="f [cycles/pix]", ylim=(-1.05, 1.05))
    for a in ax[3:]:
        a.title.set_fontsize(9)
        a.grid(alpha=0.3, which="both")
    return fig


_GRID = {1: (1, 1), 2: (1, 2), 3: (1, 3), 4: (2, 2)}  # subplot layout per number of images


def plot_single(imgs, titles, cfg, overlays=None):
    """Images in original units, each with its own color scale; zoom/pan shared if shapes match."""
    overlays = overlays or [None] * len(imgs)
    nr, nc = _GRID[len(imgs)]
    share = len({i.shape for i in imgs}) == 1
    fig = Figure(figsize=(6 * nc, 5.5 * nr), layout="constrained")
    ax0 = None
    for k, (img, title, ov) in enumerate(zip(imgs, titles, overlays)):
        ax = fig.add_subplot(nr, nc, k + 1, sharex=ax0 if share else None, sharey=ax0 if share else None)
        ax0 = ax0 or ax
        _imshow(fig, ax, img, title, cfg["cmap"], _clim(cfg["clip_percentiles"], img), overlay=ov)
        ax.set(xlabel="x [pix]", ylabel="y [pix]")
    return fig


def plot_batch(rows, method, scale=None):
    """Batch trends vs image index (0 = reference): displacement, widths, flux ratio, rms_rel."""
    u = "\u00b5m" if scale else "pix"
    k = scale or 1.0
    idx = np.array([r["index"] for r in rows])
    get = lambda key, f=1.0: np.array([r.get(key, np.nan) for r in rows], float) * f
    fig = Figure(figsize=(12, 8), layout="constrained")
    ax = [fig.add_subplot(2, 2, i + 1) for i in range(4)]
    ax[0].plot(idx, get("dpos_x", k), "o-", label="\u0394x")
    ax[0].plot(idx, get("dpos_y", k), "s-", label="\u0394y")
    ax[0].set(title=f"{'Fit center' if method == 'gauss' else 'Centroid'} displacement B \u2212 A", ylabel=f"[{u}]")
    wkeys = (("fwhm_x", "FWHM x"), ("fwhm_y", "FWHM y")) if method == "rainer" else (
        ("d_sigma_x", "d\u03c3x"), ("d_sigma_y", "d\u03c3y"), ("d_sigma", "d\u03c3"))
    for key, lab in wkeys:
        ax[1].plot(idx, get(key, k), "o-", label=lab)
    ax[1].set(title="Spot widths", ylabel=f"[{u}]")
    ax[2].plot(idx, get("flux_ratio"), "o-")
    ax[2].axhline(1, ls=":", c="k")
    ax[2].set(title="Total flux ratio B / A")
    ax[3].plot(idx, get("rms_rel"), "o-")
    ax[3].set(title="rms_rel of A \u2212 B (after registration)")
    for a in ax:
        a.set_xlabel("image index (0 = reference)")
        a.grid(alpha=0.3)
        a.title.set_fontsize(9)
        if a.get_legend_handles_labels()[0]:
            a.legend(fontsize=8)
    return fig
