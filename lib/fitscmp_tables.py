"""Side-panel tables: raw image statistics and spot metrics, with definitions for the user."""
from .fitscmp_analysis import image_stats
from .fitscmp_spot import ROUND_LIMIT, SPOT_METHODS, SpotError, spot_metrics

SECTION = "section"  # row kind for section titles

STAT_ROWS = (("min", "Min"), ("max", "Max"), ("mean", "Mean"), ("median", "Median"),
             ("std", "Std"), ("sum", "Sum"), ("n_valid", "Valid pixels"), ("n_nan", "NaN pixels"))

# (key, label, delta, is_length): delta is "diff" (B - A), "ratio" (B / A) or None;
# lengths and positions get the unit in the label and are converted to um if requested
_CENTROID = (("centroid_x", "Centroid x", "diff", True), ("centroid_y", "Centroid y", "diff", True))
_PEAK = (("peak_value", "Peak value", "diff", False), ("peak_x", "Peak position x", "diff", True),
         ("peak_y", "Peak position y", "diff", True))
_SHAPE = (("phi_deg", "Azimuth \u03c6 [deg]", "diff", False), ("ellipticity", "Ellipticity", "diff", False),
          ("round", f"Round (ellipticity > {ROUND_LIMIT})", None, False))
SPOT_ROWS = {
    "iso": (("background", "Background", "diff", False), ("background_rms", "Background rms", "diff", False))
    + _CENTROID + _PEAK
    + (("total_flux", "Total flux", "ratio", False),
       ("d_sigma_x", "Width d\u03c3x", "diff", True), ("d_sigma_y", "Width d\u03c3y", "diff", True),
       ("d_sigma", "Width d\u03c3", "diff", True))
    + _SHAPE
    + (("n_saturated", "Saturated pixels", "diff", False),
       ("area_in_frame", "Integration area in frame", None, False), ("n_iter", "Iterations", None, False)),
    "gauss": (("background", "Background (fit)", "diff", False),
              ("center_x", "Fit center x", "diff", True), ("center_x_err", "Fit center x error", None, True),
              ("center_y", "Fit center y", "diff", True), ("center_y_err", "Fit center y error", None, True))
    + _PEAK
    + (("amplitude", "Fit amplitude", "diff", False), ("total_flux", "Total flux (fit)", "ratio", False),
       ("d_sigma_x", "Width d\u03c3x (4\u03c3)", "diff", True), ("d_sigma_x_err", "Width d\u03c3x error", None, True),
       ("d_sigma_y", "Width d\u03c3y (4\u03c3)", "diff", True), ("d_sigma_y_err", "Width d\u03c3y error", None, True),
       ("d_sigma", "Width d\u03c3", "diff", True),
       ("fwhm_x", "FWHM x", "diff", True), ("fwhm_y", "FWHM y", "diff", True))
    + _SHAPE
    + (("n_saturated", "Saturated pixels (excluded)", "diff", False), ("n_fit", "Fitted pixels", None, False),
       ("residual_ratio", "Fit residual / noise", None, False), ("converged", "Converged", None, False)),
    "rainer": (("background", "Background (p10)", "diff", False),)
    + _CENTROID + _PEAK
    + (("total_flux", "Total flux (clipped)", "ratio", False),
       ("fwhm_x", "FWHM x", "diff", True), ("fwhm_y", "FWHM y", "diff", True),
       ("n_saturated", "Saturated pixels", "diff", False)),
}

_PEAK_DEF = "Peak position: position of the brightest pixel (mean position if several share the maximum). "
DEFINITIONS = {
    "iso": ("Centroid: flux-weighted mean position (first moments, ISO 11146) of the "
            "background-subtracted signal in the integration area. " + _PEAK_DEF +
            "d\u03c3x, d\u03c3y: second-moment widths (4\u03c3) along the principal axes; \u03c6: angle between "
            "the x axis and the principal axis closer to it (undefined for round beams); "
            "ellipticity: min/max width ratio. Integration area: 3 x widths."),
    "gauss": ("Fit center: center of an elliptical 2D Gaussian + constant background fitted by least squares "
              "on the ISO integration area (not the centroid). " + _PEAK_DEF +
              "d\u03c3x, d\u03c3y = 4\u03c3 and FWHM along the principal axes (ISO convention for x, y and \u03c6). "
              "Saturated pixels are excluded from the fit. Errors: 1\u03c3. Fit residual / noise >> 1 "
              "(tens) indicates a non-Gaussian profile, e.g. a near field."),
    "rainer": ("Rainer (p10) method (Acquisition Manager): background = 10th percentile of the frame, "
               "negative residuals clipped, centroid and FWHM (2.3548 \u03c3, along the frame axes) from the "
               "whole frame. Kept for comparison with existing files; biased when the frame is much larger "
               "than the spot. " + _PEAK_DEF),
}


def measure_spots(imgs, method="iso", saturation=None):
    """[(metrics or None, error message or None)] for each image.

    `saturation`: None, a level, or a list with one level or boolean mask per image.
    """
    sats = saturation if isinstance(saturation, (list, tuple)) else [saturation] * len(imgs)
    out = []
    for img, sat in zip(imgs, sats):
        try:
            out.append((spot_metrics(img, method, sat), None))
        except SpotError as e:
            out.append((None, str(e)))
    return out


def _display(m, method):
    """Copy of the metrics with the display conventions (undefined azimuth, iteration text)."""
    if m is None:
        return None
    m = dict(m)
    if m.get("round"):
        m["phi_deg"] = None  # azimuth undefined for round beams
    if method == "iso":
        m["n_iter"] = f"{m['n_iter']}" + ("" if m["converged"] else " (not converged)")
    return m


def _is_num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _delta(a, b, kind):
    """B - A, or B / A as a string, for numeric values; None otherwise."""
    if kind is None or not (_is_num(a) and _is_num(b)):
        return None
    if kind == "ratio":
        return f"\u00d7{b / a:.5g}" if a else None
    return b - a


def build_table(imgs, names, spots, method="iso", saturation=None, delta=False, scale=None, extra_stats=()):
    """Rows for the result side panel and a note with definitions and errors.

    `spots` comes from measure_spots(). Rows are (SECTION, title) or (label, values); with
    delta=True (two images) each value row has a third column B - A (B / A for fluxes).
    With a pixel scale `scale` [um/pix], lengths and positions are given in um.
    `extra_stats`: additional (label, values) rows for the statistics section.
    """
    unit = "\u00b5m" if scale else "pix"
    rows = [(SECTION, "Image statistics")]
    stats = [image_stats(i) for i in imgs]
    for key, lab in STAT_ROWS:
        rows.append((lab, [s.get(key) for s in stats] + ([None] if delta else [])))
    for lab, vals in extra_stats:
        rows.append((lab, list(vals) + ([None] if delta else [])))
    shown = [_display(m, method) for m, _ in spots]
    rows.append((SECTION, f"Spot: {SPOT_METHODS[method]}"))
    for key, lab, kind, is_len in SPOT_ROWS[method]:
        if key == "n_saturated" and saturation is None:
            continue
        vals = [m.get(key) if m else None for m in shown]
        if is_len and scale:
            vals = [v * scale if _is_num(v) else v for v in vals]
        label = f"{lab} [{unit}]" if is_len else lab
        rows.append((label, vals + ([_delta(vals[0], vals[1], kind)] if delta else [])))
    errors = [f"{name}: spot not measured ({err})." for name, (_, err) in zip(names, spots) if err]
    coords = ("All values in original units; positions are 0-based pixel indices (x = column, y = row; "
              "add 1 for FITS/DS9)" + (f" multiplied by the pixel scale {scale:g} \u00b5m/pix." if scale else "."))
    return rows, " ".join([DEFINITIONS[method], coords] + errors)


# ---------------------------------------------------------------- batch table
def batch_columns(method, saturation=None, hot_filter=False, scale=None):
    """Batch columns (key, display label, export name, is_length) for the selected spot method."""
    u, eu = ("\u00b5m", "um") if scale else ("pix", "pix")
    pos = ("center_x", "center_y") if method == "gauss" else ("centroid_x", "centroid_y")
    pname = "Fit center" if method == "gauss" else "Centroid"
    cols = [("index", "#", "index", False), ("image", "Image", "image", False),
            (pos[0], f"{pname} x [{u}]", f"{pos[0]}_{eu}", True), (pos[1], f"{pname} y [{u}]", f"{pos[1]}_{eu}", True),
            ("dpos_x", f"\u0394x B\u2212A [{u}]", f"dpos_x_{eu}", True), ("dpos_y", f"\u0394y B\u2212A [{u}]", f"dpos_y_{eu}", True)]
    if method in ("iso", "gauss"):
        cols += [("d_sigma_x", f"d\u03c3x [{u}]", f"d_sigma_x_{eu}", True), ("d_sigma_y", f"d\u03c3y [{u}]", f"d_sigma_y_{eu}", True),
                 ("d_sigma", f"d\u03c3 [{u}]", f"d_sigma_{eu}", True), ("phi_deg", "\u03c6 [deg]", "phi_deg", False),
                 ("ellipticity", "Ellipticity", "ellipticity", False)]
    if method in ("gauss", "rainer"):
        cols += [("fwhm_x", f"FWHM x [{u}]", f"fwhm_x_{eu}", True), ("fwhm_y", f"FWHM y [{u}]", f"fwhm_y_{eu}", True)]
    cols += [("total_flux", "Flux", "total_flux", False), ("flux_ratio", "Flux B/A", "flux_ratio", False),
             ("peak_value", "Peak", "peak_value", False)]
    if method == "gauss":
        cols += [("residual_ratio", "Fit resid./noise", "fit_residual_ratio", False)]
    if saturation is not None:
        cols += [("n_saturated", "Saturated", "n_saturated", False)]
    if hot_filter:
        cols += [("n_hot", "Hot px", "n_hot", False)]
    cols += [("shift_dy", "Shift dy [pix]", "shift_dy_pix", False), ("shift_dx", "Shift dx [pix]", "shift_dx_pix", False),
             ("rms_rel", "rms_rel", "rms_rel", False), ("pearson", "Pearson", "pearson", False),
             ("frc_f50", "FRC 0.5 [cyc/pix]", "frc_f50_cycles_per_pix", False), ("error", "Error", "error", False)]
    return cols


def batch_values(rows, cols, scale=None):
    """Table values (one list per image), lengths converted to um if a scale is given."""
    out = []
    for r in rows:
        vals = []
        for key, _, _, is_len in cols:
            v = r.get(key)
            if key == "error":
                v = "; ".join(e for e in (r.get("error", ""), r.get("spot_error", "")) if e) or ("reference" if r.get("reference") else "")
            elif is_len and scale and _is_num(v):
                v = v * scale
            vals.append(v)
        out.append(vals)
    return out
