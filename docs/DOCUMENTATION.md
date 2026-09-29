# fitscmp — Technical Documentation

Tool for the quantitative comparison of two 2D images stored in one or more FITS files, both in real space (direct subtraction) and in Fourier space (power spectra, radial profiles, Fourier Ring Correlation).

---

## Contents

1. [Overview](#1-overview)
2. [Project structure](#2-project-structure)
3. [Requirements and installation](#3-requirements-and-installation)
4. [Quick start](#4-quick-start)
5. [Configuration (`fitscmp.yaml`)](#5-configuration-fitscmpyaml)
6. [Graphical interface](#6-graphical-interface)
7. [FITS I/O (`lib/fitscmp_io.py`)](#7-fits-io-libfitscmp_iopy)
8. [Analysis methods (`lib/`)](#8-analysis-methods-lib)
9. [Output figures](#9-output-figures)
10. [Interpretation notes](#10-interpretation-notes)
11. [API reference](#11-api-reference)
12. [Synthetic test data (`make_test_fits.py`)](#12-synthetic-test-data-make_test_fitspy)
13. [Known limitations](#13-known-limitations)
14. [References](#14-references)

---

## 1. Overview

The processing chain for a comparison is:

```
FITS files ──► list_images() ──► pair selection (A, B) ──► load_image() ──► calibrate() ──► float64 arrays
                                                                                               │
                                    [estimate_shift()] ──► apply_shift(B)  ◄───────────────────┘
                                                                 │
                              normalize(A), normalize(B)  ◄──────┘
                                         │
                ┌────────────────────────┴────────────────────────┐
         difference()                              fourier_compare() + fitscmp_mask
   A−B map + scalar statistics            2D log power spectra, log ratio, radial profiles, FRC
                │                                                  │
        plot_subtraction()                                   plot_fourier()
                └──────────────► ResultWindow (two tabs) ◄─────────┘
```

Image B is always compared against image A (A is the reference in asymmetric quantities such as `rms_rel` and the ratio `P(A)/P(B)`). Image statistics and spot metrics (side table, [§8.5](#85-spot-analysis-libfitscmp_spotpy)) are computed on the images as loaded, in original units; the single-image view ([§6.2](#62-image-selection)) shows 1 to 4 images without comparison, and the batch ([§6.6](#66-batch-comparison)) compares a series of images with a reference.

The processing is implemented in `lib/fitscmp_pipeline.py` (`run_single`, `run_comparison`, `run_batch`), independent of the GUI: the GUI collects the options into a `Settings` object, runs the pipeline in the background ([§6.5](#65-background-execution)) and displays the results. The same functions can be used from scripts ([§11](#11-api-reference)).

---

## 2. Project structure

```
Compare_fits/
├── fitscmp_gui.py          # Tkinter GUI, entry point
├── fitscmp.yaml            # Default settings
├── make_test_fits.py       # Generator of synthetic test FITS files
├── requirements.txt        # Python dependencies (pip)
├── CHANGELOG.md            # History of changes
├── .gitignore
├── lib/
│   ├── __init__.py
│   ├── fitscmp_io.py       # FITS discovery, loading, header extraction
│   ├── fitscmp_analysis.py # Normalization, subtraction statistics, Fourier analysis
│   ├── fitscmp_mask.py     # Masked pixels in the FFT: apodized weights, inpainting
│   ├── fitscmp_registration.py # Sub-pixel shift estimate and Fourier-shift resampling
│   ├── fitscmp_spot.py     # Laser-spot metrics: ISO 11146, Rainer (p10)
│   ├── fitscmp_fit.py      # Elliptical Gaussian fit of the spot
│   ├── fitscmp_calib.py    # Master dark subtraction, hot-pixel filter
│   ├── fitscmp_pipeline.py # Processing pipeline: single view, comparison, batch (no GUI)
│   ├── fitscmp_export.py   # Export of results: CSV / Parquet tables, JSON metadata
│   ├── fitscmp_tables.py   # Side-panel tables and definitions
│   ├── fitscmp_plots.py    # Matplotlib figures (backend-independent)
│   └── fitscmp_windows.py  # Tk windows: dialogs, header browser, result window
├── tests/
│   ├── test_spot.py        # Validation of the spot metrics on simulated spots
│   ├── test_calib.py       # Validation of dark subtraction and hot-pixel filter
│   └── test_pipeline.py    # Validation of batch, process pool, cancel and export
└── docs/
    └── DOCUMENTATION.md    # This file
```

All modules in `lib/` except `fitscmp_windows.py` are independent of Tkinter (the plotting module uses only `matplotlib.figure`) and can be used from scripts or notebooks (see [§11](#11-api-reference)).

---

## 3. Requirements and installation

| Package | Use | Notes |
|---|---|---|
| Python 3 | — | Developed with Python 3.13 |
| NumPy | arrays, FFT | |
| Astropy | FITS I/O | `astropy.io.fits` |
| Matplotlib | figures | ≥ 3.5 (uses `Figure(layout="constrained")`) and the `TkAgg` backend |
| SciPy | distance transform, mask shifting | |
| scikit-image | shift estimate (`phase_cross_correlation`) | |
| PyYAML | configuration | |
| pandas, pyarrow | Parquet export (optional) | only for `export_format: parquet` |
| Tkinter | GUI | Part of the standard library, but on some Linux distributions it is a separate system package (e.g. `python3-tk`) |

```bash
pip install -r requirements.txt
python tests/test_spot.py     # or: python -m pytest tests
python tests/test_calib.py
python tests/test_pipeline.py
```

Verified with Python 3.12, NumPy 2.4, Astropy 8.0, SciPy 1.17, scikit-image 0.26, Matplotlib 3.10 and Tk 8.6 (GUI tested under a virtual X server).

---

## 4. Quick start

```bash
cd Compare_fits
python make_test_fits.py      # writes test_single.fits, test_mef.fits, test_spot.fits and test_dark.fits
python fitscmp_gui.py
```

In the GUI: **Add files...** → select the two test files → **Compare** → choose images A and B in the dialog → inspect the *Subtraction* and *Fourier* tabs. To look at images without comparing them, check *Single images* (the button becomes **Show**) or add a file containing a single image.

The GUI imports `lib` as a package relative to the working directory, so it must be launched from the project root (or with the root on `PYTHONPATH`). The configuration file, instead, is always searched next to `fitscmp_gui.py`.

---

## 5. Configuration (`fitscmp.yaml`)

The YAML file is loaded on top of the built-in defaults (`DEFAULTS` in `fitscmp_gui.py`); missing keys keep the default value, and a missing file is not an error.

| Key | Default | GUI control | Description |
|---|---|---|---|
| `normalization` | `none` | yes | Flux normalization: `none`, `mean`, `median`, `max`, `sum` ([§8.1](#81-normalization)). Invalid values fall back to `none`. |
| `window` | `true` | yes | Apply a 2D Hann window before the FFT. |
| `radial_bins` | `64` | yes (8–1024) | Number of annuli between 0 and 0.5 cycles/pixel. |
| `clip_percentiles` | `[1.0, 99.0]` | no | Percentiles used for the color limits of all maps. |
| `cmap` | `viridis` | no | Colormap for images and log power spectra. |
| `diff_cmap` | `RdBu_r` | no | Colormap for signed maps (A−B, log ratio). |
| `nan_mode` | `interp` | yes | Pixels not valid in both images, Fourier analysis only: `interp` or `mask` ([§8.3](#83-fourier-analysis)). Invalid values fall back to `interp`. |
| `mask_taper` | `16` | yes (0–64) | Width [pix] of the cosine apodization around masked regions; 0 disables it (negative values act as 0). |
| `interp_sigma` | `1.0` | no | Standard deviation [pix] of the Gaussian inpainting kernel (`interp` mode). |
| `registration_upsample` | `100` | no | Upsampling factor $u$ of the shift estimate (precision $1/u$ pix). |
| `spot_method` | `iso` | yes | Spot analysis method: `iso` (ISO 11146), `gauss` (Gaussian fit) or `rainer` (Rainer p10) ([§8.5](#85-spot-analysis-libfitscmp_spotpy)). Invalid values fall back to `iso`. |
| `pixel_scale` | `null` | *Units µm* | Detector pixel size [µm/pix]; enables the *Units µm* option (disabled if `null`). |
| `dark_frame` | `null` | *Master dark...* | Master dark FITS file loaded at start (first 2D image), subtracted from all images ([§8.6](#86-frame-calibration-libfitscmp_calibpy)). |
| `hot_pixel_nsigma` | `5` | shown | Detection threshold of the hot-pixel filter, in units of the local noise. |
| `saturation_level` | `null` | shown | Pixel value counted as saturated (e.g. 4095 for 12-bit data); `null` disables the check. |
| `export_format` | `csv` | no | Format of the exported tables: `csv` or `parquet` (Brotli-compressed; needs pandas and pyarrow) ([§9.5](#95-exported-files)). |
| `batch_workers` | `auto` | no | Processes used by the batch: `auto` = min(4, CPUs − 1), or an integer; 1 = no process pool ([§6.5](#65-background-execution)). |
| `expand_cubes` | `true` | no | Treat each plane of a 3D HDU as a separate 2D image. |
| `start_dir` | `"."` | no | Initial directory of the file dialog (relative paths refer to the current working directory). |

The values in the GUI are initialized from the configuration and can be changed per comparison; they are not written back to the file. Apart from `normalization`, `nan_mode` and `spot_method`, values are not validated.

---

## 6. Graphical interface

### 6.1 Main window (`App`)

- **File list**: FITS files added by the user (duplicates are ignored). Extended selection is supported.
- **Add files...**: file dialog filtered on `*.fits *.fit *.fts *.fz *.gz`.
- **Remove / Clear**: remove the selected files / all files.
- **Inspect**: opens the header browser for the selected files.
- **Compare** / **Show**: collects all 2D images from the selected files and starts the comparison, or the single-image view if *Single images* is checked or only one image is available ([§6.2](#62-image-selection)). The button label follows the mode.
- **Batch...**: compares all images of the selected files with a reference chosen in a dialog ([§6.6](#66-batch-comparison)).
- **Status line** (bottom): messages, progress bar and **Cancel** for the running job ([§6.5](#65-background-execution)).
- **Options**: *Single images* (view mode), normalization, Hann window, number of radial bins, *Swap A/B* ([§6.2](#62-image-selection)). Normalization, FFT and registration settings apply to the comparison only.
- **Masked pixels (FFT)**: NaN mode and taper width ([§8.3](#83-fourier-analysis)).
- **Calibration** (applied to all images, both views, [§8.6](#86-frame-calibration-libfitscmp_calibpy)): *Subtract dark* (enabled once a master dark is loaded), *Master dark...* (select the file), *Clear*, the label of the loaded dark, and *Hot-pixel filter*.
- **Spot analysis**: method (*ISO 11146*, *Gaussian fit* or *Rainer (p10)*, [§8.5](#85-spot-analysis-libfitscmp_spotpy)); *Overlay* (spot markers on the images, [§9.4](#94-spot-overlay)); *Units µm* (positions and widths in the table converted with `pixel_scale`; the image axes stay in pixels, consistent with the overlay and the cursor readout; disabled if `pixel_scale` is not set); the saturation level from the configuration is shown. Applies to both views.
- **Registration**: shift (dy, dx) applied to B, *Auto (estimate per comparison)* (estimate the shift for each comparison, write it in the fields with 3 decimals and apply it; in the batch, estimated for each image), *Reset* (set the shift to 0). To refine an automatic estimate by hand, uncheck *Auto*, edit the fields and compare again. With *Swap A/B* the shift is still applied to the image that ends up as B ([§8.4](#84-registration-libfitscmp_registrationpy)).

*Inspect* and *Compare*/*Show* operate on the selected files, or on all files if nothing is selected.

### 6.2 Image selection

All 2D images found in the target files are enumerated (see [§7](#7-fits-io-libfitscmp_iopy)).

- No image → error.
- Exactly 1 image → single-image view, whatever the state of *Single images*.

**Single-image view** (*Single images* checked): if more than one image is available, a modal `SelectDialog` lets the user choose 1 to 4 of them (Ctrl/Shift for multiple selection).

**Comparison** (*Single images* unchecked, at least 2 images):

- Exactly 2 images → they are used directly, in list order (first = A), or in reverse order if *Swap A/B* is checked.
- More than 2 → a modal `PairDialog` lets the user choose A and B; they must differ. The default selection is #0 → A, #1 → B, inverted if *Swap A/B* is checked; the user's choice in the dialog always prevails. Images are identified by the label `file [HDU n 'EXTNAME', plane k] NYxNX`.

The two images must have the same shape; otherwise an error is shown after loading.

### 6.3 Header browser (`HeaderWindow`)

One list entry per (file, HDU); selecting it shows the full header as text. Files that cannot be opened appear as an `ERROR` entry with the exception message.

### 6.4 Result window (`ResultWindow`)

Shows a header with the image labels, the figure tabs, each with a Matplotlib navigation toolbar (zoom, pan, save figure), and a side table in **original units**, with two sections: *Image statistics* (min, max, mean, median, standard deviation, sum, number of valid and NaN pixels, on the finite pixels) and *Spot* (metrics of the selected method, [§8.5](#85-spot-analysis-libfitscmp_spotpy)). A note under the table defines centroid, peak position and widths, recalls the coordinate convention and reports images whose spot could not be measured. The header also reports the calibration applied.

**Pixel readout.** Above the tabs, a line shows, for the pixel under the cursor on any image panel, its 0-based coordinates (also in µm with *Units µm*) and the value of **every** image at that pixel: `#k` in the single-image view (images of the same shape as the hovered one), and A, B (as displayed, `B(shifted)` if a shift is applied) and A−B (normalized, as in the panel) in the comparison. A left click pins the current readout (up to 5 pins, *Clear pins* removes them); clicks are ignored while zoom or pan are active. The readout is off on the Fourier spectra. Values are the calibrated images in original units.

- Comparison: header with A, B and the shift applied to B (marked `[auto]` if estimated); tabs *Subtraction* and *Fourier*; table columns A, B and B − A. The table always refers to the images **as recorded** (B before the shift). B − A is the difference, except for the total flux (ratio B/A, shown as ×value); it is empty for counts and flags. Comparing the centroid difference with the *Auto* shift gives two independent estimates of the displacement: they differ when the spot shape changes, since the cross-correlation aligns the whole profile.
- Single-image view: header with the labels `#k`; tab *Images*; one table column per image.

**Export...** (top right) writes the numerical results ([§9.5](#95-exported-files)). Multiple result windows can be open at the same time.

When the window is mapped, each figure is resized to its widget. This works around a Matplotlib TkAgg ordering issue on Linux screens whose DPI differs from 96 (device pixel ratio ≠ 1), where the figure would otherwise stay a few percent larger than a canvas that cannot grow, and its right side would be clipped.

### 6.5 Background execution

Computations run outside the GUI thread, so that the interface stays responsive; Tk is used only from the GUI thread, which polls the worker through a queue.

- **Compare / Show** run in a worker thread. **Batch** runs in a worker thread that distributes the images to a process pool (`concurrent.futures.ProcessPoolExecutor`, *spawn* start method, which is safe with a Tk application); `batch_workers` = `auto` uses min(4, CPUs − 1) processes, an integer sets their number, and 1 processes the images sequentially in the worker thread. The process pool gives identical results to the sequential run (`tests/test_pipeline.py`).
- While a job runs, *Compare/Show* and *Batch...* are disabled (one job at a time), the progress bar at the bottom of the main window runs (determinate for the batch, with the image count in the status line) and **Cancel** is enabled.
- **Cancel** stops a batch between two images (in the process pool, images already being processed are completed, the others are dropped) and no window is opened. For a single comparison or view, the computation cannot be interrupted: its result is discarded when it ends.
- The options are read when the job starts (`Settings`, [§11](#11-api-reference)): changing them during a job affects the next one only.

### 6.6 Batch comparison

**Batch...** takes all 2D images of the target files; a dialog (`BatchDialog`) selects the **reference** (A, index 0), and every other image (B, indices 1…N in list order) is compared with it, using the current options (calibration, spot method, normalization, Fourier settings, registration: with *Auto* the shift is estimated for each image, otherwise the fixed shift of the fields is applied to all). Images with a shape different from the reference get their spot metrics only, with an error message.

The result window (`BatchWindow`) has a **Table** tab, one row per image:

| Columns | Content |
|---|---|
| #, Image | index and label |
| Centroid (or Fit center) x, y | spot position as recorded |
| Δx, Δy B−A | position minus the reference position, **before** registration |
| dσx, dσy, dσ, φ, Ellipticity / FWHM x, y | widths and shape of the selected method; φ is blank for round beams |
| Flux, Flux B/A, Peak | total flux, ratio to the reference, peak value |
| Fit resid./noise, Saturated, Hot px | only for the Gaussian fit, with a saturation level, with the hot-pixel filter |
| Shift dy, dx [pix] | registration applied before the difference |
| rms_rel, Pearson | statistics of A − B ([§8.2](#82-real-space-difference)), after registration and normalization |
| FRC 0.5 [cyc/pix] | first spatial frequency where the Fourier Ring Correlation with the reference drops below 0.5 (linear interpolation between annuli; the first annulus center if already below; blank if never) |
| Error | reason why a quantity is missing |

and a **Trends** tab with four panels versus the image index: displacement Δx, Δy; widths; flux ratio; rms_rel. With *Units µm*, positions and widths are in µm. A note under the table recalls the definitions.

---

## 7. FITS I/O (`lib/fitscmp_io.py`)

### 7.1 Image discovery

`list_images(path, expand_cubes=True)` scans all HDUs of a file and returns a list of `ImageRef`. Only header information is read; no pixel data is loaded.

| HDU type | Dimensionality | Result |
|---|---|---|
| `PrimaryHDU`, `ImageHDU`, `CompImageHDU` | 2D | one image |
| same | 3D, `expand_cubes=True` | one image per plane along the first NumPy axis (FITS `NAXIS3`) |
| same | 3D, `expand_cubes=False` | ignored |
| same | empty, 1D, ≥ 4D | ignored |
| tables and other HDUs | — | ignored |

Tile-compressed (`.fz`) and gzip-compressed files are handled by Astropy.

### 7.2 `ImageRef`

Immutable dataclass pointing to one 2D image:

| Field | Type | Meaning |
|---|---|---|
| `path` | `str` | file path |
| `hdu` | `int` | HDU index |
| `plane` | `int \| None` | plane index for cubes, `None` for 2D HDUs |
| `extname` | `str` | `EXTNAME` (empty for unnamed and `PRIMARY` HDUs) |
| `shape` | `tuple` | `(ny, nx)` of the 2D image |

`label()` returns the human-readable identifier used in the GUI.

### 7.3 Loading

`load_image(ref)` reopens the file and returns the image as a native-endian `float64` copy. Astropy applies `BSCALE`/`BZERO` (the file is opened with the default `memmap=None`, since an explicit `memmap=True` makes Astropy refuse scaled data). The copy decouples the array from the file, which is closed on return.

### 7.4 Headers

`read_headers(path)` returns `[(title, header_text), ...]` for every HDU, with `title = "HDU i: NAME (HDUClass)"`.

---

## 8. Analysis methods (`lib/`)

Notation: $A$, $B$ are the two images, $\Omega$ the set of pixels that are finite in both, $N_\Omega = |\Omega|$.

### 8.1 Normalization

Each image is divided by its own scalar statistic, computed on the **common mask** $\Omega$ (pixels finite in both images), so that both factors refer to the same set of pixels:

$$A' = A / s(A_\Omega), \qquad B' = B / s(B_\Omega), \qquad s \in \{\text{mean},\ \text{median},\ \text{max},\ \text{sum}\}$$

A non-finite or zero factor raises `ValueError`. With `none`, images are unchanged. All following quantities (statistics, spectra, FRC) are computed on the normalized images.

Notes: since both sums run over the same $N_\Omega$ pixels, `mean` and `sum` differ only by the common factor $N_\Omega$, and all scale-free results (`rms_rel`, `pearson`, spectral ratios, FRC) are identical; `max` is sensitive to hot pixels and cosmic rays. Pixels finite in only one image are normalized but do not contribute to the factor. `normalize()` called without a mask (API use) computes the statistic on all finite pixels of the image.

### 8.2 Real-space difference

$D = A - B$ (NaN where either input is NaN). Statistics over $\Omega$:

| Key | Definition |
|---|---|
| `mean` | $\langle D \rangle$ |
| `std` | $\sigma_D$ (population, `ddof=0`) |
| `rms` | $\sqrt{\langle D^2 \rangle}$ |
| `max_abs` | $\max \lvert D \rvert$ |
| `rms_rel` | $\sqrt{\langle D^2 \rangle} \,/\, \sigma_A$ — relative to the standard deviation of A (independent of any pedestal) |
| `pearson` | Pearson correlation coefficient of $A$ and $B$ |
| `n_valid` | $N_\Omega$ |
| `masked` | $1 - N_\Omega / (N_x N_y)$, fraction of pixels not finite in both images |

For constant images `rms_rel` and `pearson` are undefined and returned as NaN/inf (no warning is emitted).

### 8.3 Fourier analysis

#### Masked pixels and weight map (`lib/fitscmp_mask.py`)

A, B and $D = A - B$ are transformed with **the same weight map**

$$W = H \cdot T \cdot M,$$

so that masking affects the three spectra in the same way:

- $M$ is the common validity mask (pixels finite in both images), computed after the optional inpainting:
  - `nan_mode = "interp"` (default): pixels invalid in A **or** B are set to NaN in both images and filled by normalized convolution with a Gaussian kernel of standard deviation `interp_sigma` (`astropy.convolution.interpolate_replace_nans`, kernel size $8\sigma + 1$). Holes larger than the kernel support stay invalid and are masked. Inpainting is used only for the Fourier analysis; real-space statistics always use $\Omega$.
  - `nan_mode = "mask"`: no filling.
- $T$ is a cosine apodization of the mask edges: $T = \tfrac{1}{2}\left[1 - \cos(\pi d/\tau)\right]$ for $d < \tau$ and $T = 1$ otherwise, where $d$ is the Euclidean distance to the nearest invalid pixel (`scipy.ndimage.distance_transform_edt`) and $\tau$ = `mask_taper`. With $\tau = 0$ the mask is binary. $T$ does not taper the image borders, which are handled by $H$.
- $H$ is the separable 2D Hann window $h_i h_j$ (`np.hanning`, symmetric, zero at the edges) [Harris 1978], or 1 if the window is disabled.

#### Preprocessing (`_centered_fft`)

For each input image $x \in \{A, B, D\}$:

1. subtract the mean computed on the pixels with $W > 0$;
2. set the pixels with $W = 0$ to 0 and multiply by $W$;
3. compute the 2D FFT and shift the zero frequency to the center.

#### Power spectrum normalization

$$P(\mathbf{k}) = \frac{\lvert \hat{x}(\mathbf{k}) \rvert^2}{\sum_{ij} W_{ij}^2}$$

which reduces to $N_x N_y$ without window and mask. By Parseval's theorem, the mean of $P$ over the whole frequency plane equals the weighted variance $\sum (x W)^2 / \sum W^2$, so spectra obtained with different windows and masks are on the same scale.

#### Pseudo-spectra and choice of defaults

With a mask, $P$ is a *pseudo-spectrum*: the true spectrum convolved with $\lvert \hat{W} \rvert^2$ (mode coupling), which is not deconvolved [Hivon et al. 2002]. Since $W$ is the same for A and B, the coupling is identical in both, and ratios and FRC are preserved to first order. Sharp mask edges, however, leak power into all frequencies *identically* in A and B; this leakage is correlated and inflates the FRC. Isolated bad pixels are the worst case, because each one is a sharp hole. This is why the default fills small defects by inpainting and apodizes the residual holes.

The defaults were chosen on simulations: a speckle pattern with a gradient, independent noise in A and B, B blurred (Gaussian, $\sigma = 1$ pix), and NaN defects in B only (disk of radius 20 pix, one bad column, 0.5% random pixels). The table gives the median over 15 realizations of the maximum deviation from the defect-free case:

| Method | max $\lvert \Delta \log_{10} (\bar{P}_A/\bar{P}_B) \rvert$ | max $\lvert \Delta\,\mathrm{FRC} \rvert$ |
|---|---|---|
| Original (independent mean fill) | 0.27 | 0.23 |
| `mask`, taper 16 | 0.14 | 0.21 |
| `interp`, taper 8 | 0.06 | 0.09 |
| `interp`, taper 16 (**default**) | 0.05 | 0.08 |

With the large hole alone, `mask` and `interp` are equivalent (≈ 0.04 / 0.05 with taper 8). On `test_mef.fits` (steeper spectrum, see [§12](#12-synthetic-test-data-make_test_fitspy)), both new modes are exact below $f \approx 0.15$ cycles/pixel, while the original method deviates by ≈ 0.03; a residual leakage remains in the noise-dominated band, decreasing with the taper width.

Subtracting the $W$-weighted mean (which zeroes the DC term of the weighted image) was also tested and rejected: on white noise it biases the first radial bin to 0.85 instead of 0.96 (expected value 1).

#### Frequency grid

Frequencies are in cycles/pixel from `np.fft.fftfreq`, in the range $[-0.5, 0.5)$; axis 0 is $f_y$, axis 1 is $f_x$, and $f = \sqrt{f_x^2 + f_y^2}$. For non-square images the frequency sampling differs along the two axes ($\Delta f_x = 1/N_x$, $\Delta f_y = 1/N_y$).

#### 2D maps

- `logps_a`, `logps_b`: $\log_{10} P_A$, $\log_{10} P_B$;
- `log_ratio`: $\log_{10}(P_A/P_B)$.

Non-positive values become NaN.

#### Radial profiles

Bin edges $e_k = 0.5\,k/N_\text{bins}$, $k = 0 \dots N_\text{bins}$. A frequency pixel belongs to bin $k$ if $e_k \le f < e_{k+1}$; the DC pixel ($f = 0$) and all pixels with $f \ge 0.5$ (corners of the plane) are excluded. The profile is the mean of $P$ in each annulus; empty annuli give NaN. Bin centers are returned as `f`.

Profiles are computed for $P_A$, $P_B$ and $P_{A-B}$; the GUI also plots the ratio of the radial profiles $\bar{P}_A / \bar{P}_B$.

#### Fourier Ring Correlation

$$\mathrm{FRC}(k) = \frac{\sum_{\mathbf{k} \in R_k} \mathrm{Re}\!\left[\hat{A}(\mathbf{k})\,\hat{B}^*(\mathbf{k})\right]}
{\sqrt{\sum_{\mathbf{k} \in R_k} \lvert \hat{A} \rvert^2 \;\sum_{\mathbf{k} \in R_k} \lvert \hat{B} \rvert^2}}$$

where $R_k$ is the $k$-th annulus [Saxton & Baumeister 1982; van Heel & Schatz 2005]. Using the real part is exact for real images, since the imaginary contributions of $\mathbf{k}$ and $-\mathbf{k}$ cancel. FRC is 1 for identical structure at that spatial frequency, ~0 for uncorrelated content, and negative for anti-correlated content.

### 8.4 Registration (`lib/fitscmp_registration.py`)

Only translations are handled. The shift $(\Delta y, \Delta x)$ is always applied to **B**, to bring it onto A, before normalization and all other computations. Positive values move the content of B towards larger row/column indices (up and right in the displayed images).

#### Estimate (`estimate_shift`)

1. Pixels invalid in A or B are inpainted in both, as in `interp` mode ([§8.3](#83-fourier-analysis)).
2. Both images are mean-subtracted and multiplied by the same weight map $W$ (Hann window × apodized common mask, taper from the GUI).
3. The shift is the peak of the cross-correlation, refined to $1/u$ pixel by the upsampled-DFT method of Guizar-Sicairos et al. (2008), via `skimage.registration.phase_cross_correlation` with `upsample_factor` $u$ = `registration_upsample` and `normalization=None`.

Two choices were made on simulations (speckle patterns, random shifts in $[-3, 3]$ pix, 25 realizations, RMS error of the estimate):

| Setting | Band-limited signal (pupil 0.05 cycles/pix) | Broad-band signal (pupil 0.15) |
|---|---|---|
| Phase correlation (skimage default) | 0.07 – 0.25 pix | ≤ 0.006 pix |
| Plain cross-correlation (**used**) | ≤ 0.008 pix | ≤ 0.004 pix |

| Non-periodic images, plain cross-correlation | Hann window (**used**) | no window |
|---|---|---|
| pupil 0.05, no gradient / with gradient | 0.005 / 0.005 pix | 0.043 / 0.098 pix |

Phase correlation whitens the cross-power spectrum, giving equal weight to noise-dominated frequencies. Without the window, the border discontinuities of non-periodic images bias the estimate towards zero shift. The masked mode of `phase_cross_correlation` is not used because it disables sub-pixel refinement.

#### Resampling (`apply_shift`)

B is shifted with the Fourier shift theorem, i.e. multiplied by $e^{-2\pi i (f_y \Delta y + f_x \Delta x)}$ in Fourier space. Unlike spline interpolation, this does not low-pass filter the image, so the Fourier comparison is not biased at high frequency. To avoid wrap-around, the image is mirror-padded by $\lceil \lvert \Delta \rvert \rceil + 16$ pixels before the FFT and cropped afterwards. Then:

- the strips of width $\lceil \lvert \Delta y \rvert \rceil$ rows and $\lceil \lvert \Delta x \rvert \rceil$ columns entering from outside the field are set to NaN;
- NaN pixels of B are inpainted (small holes) or set to the mean (large holes) before the FFT; after the shift, every pixel that receives a contribution from an invalid pixel under linear interpolation of the mask is set to NaN (holes grow by at most one pixel).

The resulting NaN pixels are handled by the common mask like any other invalid pixel. Integer shifts reproduce `np.roll` in the interior.

### 8.5 Spot analysis (`lib/fitscmp_spot.py`)

Metrics of a single laser spot (e.g. the near or far field of a fiber output), computed on each image **in original units**. In the comparison, B is measured **before** the registration shift, so positions refer to the detector. Three methods are available (`spot_method`, GUI *Spot analysis → Method*): ISO 11146 (default), Gaussian fit, and Rainer (p10). Images are calibrated first ([§8.6](#86-frame-calibration-libfitscmp_calibpy)); the saturation check uses the raw values.

#### Definitions

| Quantity | Definition |
|---|---|
| **Centroid** $(\bar{x}, \bar{y})$ | Flux-weighted mean position (first moments) of the background-subtracted signal $E$ in the integration area: $\bar{x} = \sum E\,x / \sum E$ |
| **Fit center** $(x_0, y_0)$ | Center of the fitted Gaussian model (Gaussian fit method only); it coincides with the centroid for a noiseless Gaussian, but is a different estimator |
| **Peak position** | Position of the brightest pixel; if several pixels share the maximum (e.g. saturation), their mean position. It is not a fit and is noisy by construction: on a flat-topped near field it is an arbitrary pixel of the plateau, and the centroid is the relevant position |
| **Peak value** | Value of the brightest pixel, in original units, background **not** subtracted (to compare with the saturation level) |
| **Total flux** | $\sum E$ in the integration area |
| **Saturated pixels** | Number of pixels $\ge$ `saturation_level` (row shown only if the level is set) |

#### ISO 11146 method (default)

Beam widths, azimuth and ellipticity follow the definitions of ISO 11146-1. With the second central moments

$$\sigma_x^2 = \frac{\sum E (x-\bar{x})^2}{\sum E}, \quad \sigma_y^2 = \frac{\sum E (y-\bar{y})^2}{\sum E}, \quad \sigma_{xy}^2 = \frac{\sum E (x-\bar{x})(y-\bar{y})}{\sum E},$$

the widths along the principal axes are

$$d_{\sigma x,y} = 2\sqrt{2}\,\Bigl\{\sigma_x^2 + \sigma_y^2 \pm \gamma \bigl[(\sigma_x^2 - \sigma_y^2)^2 + 4 (\sigma_{xy}^2)^2\bigr]^{1/2}\Bigr\}^{1/2}, \qquad \gamma = \operatorname{sgn}(\sigma_x^2 - \sigma_y^2),$$

the azimuth is $\varphi = \tfrac{1}{2} \arctan\!\left[ 2\sigma_{xy}^2 / (\sigma_x^2 - \sigma_y^2) \right]$, the angle between the x axis and the principal axis closer to it ($\lvert \varphi \rvert \le 45°$), and $d_{\sigma x}$ is the width along that axis. For a round beam the width is $d_\sigma = 2\sqrt{2}\,(\sigma_x^2 + \sigma_y^2)^{1/2}$. The **ellipticity** is $\varepsilon = \min(d_{\sigma x}, d_{\sigma y}) / \max(d_{\sigma x}, d_{\sigma y})$; the beam is considered **round** if $\varepsilon > 0.87$, and then $\varphi$ is reported as undefined (—). Useful references: $d_\sigma = 4\sigma$ for a Gaussian, and $d_\sigma = 2R$ (the diameter) for a uniform disk of radius $R$.

Background and integration area follow the approach of ISO 11146-3 (background estimated outside a limited integration area, iterated), with these implementation choices:

1. **Start.** Background = sigma-clipped median (3σ) of the whole frame, with its rms. Initial spot region = the connected component with the largest flux among the pixels whose 3×3 median-filtered value exceeds background + 3 rms (isolated hot pixels do not pass the filter).
2. **Iteration.** Moments of $E$ = image − background inside the current area; new integration area = rectangle centered on the centroid, aligned with the principal axes, with sides $3\,d_{\sigma x} \times 3\,d_{\sigma y}$ (for round beams, a square of side $3\,d_\sigma$ aligned with the frame, because the azimuth is undefined); new background = sigma-clipped mean (3σ) of the pixels outside the area (kept unchanged if fewer than 100 pixels are outside). Negative residuals are **not** clipped.
3. **Convergence.** The area is a discrete set of pixels: the iteration stops when that set is stationary or repeats (limit cycle), which avoids an arbitrary tolerance; at most 50 iterations. The table reports the number of iterations and whether the area lies within the frame, as required by the standard. NaN pixels are ignored.

Why not a fixed tolerance, and why a detection-based start: on a 1024×1280 frame, $\sum r^2 \approx 3\times10^{11}$ pix², so starting from whole-frame moments, a background error of 0.01 ADU changes $\sigma^2$ by ~60%; and for round beams $\varphi$ is dominated by noise, so a rotated area changes by thousands of pixels at each iteration.

#### Rainer (p10) method

The method implemented by Monica Rainer in the Acquisition Manager camera view (`CameraView.compute_centroid`), reproduced exactly (bit-identical results, see `tests/test_spot.py`), to compare with existing results such as the `CENT_X`, `CENT_Y`, `FWHM_X`, `FWHM_Y` header keywords written by that program:

1. background = 10th percentile of all pixels of the frame;
2. image − background, negative values clipped to 0;
3. centroid and $\sigma_x$, $\sigma_y$ from the x and y projections of the **whole frame**; FWHM $= 2.3548\,\sigma$ (Gaussian profile).

On pure Gaussian noise the 10th percentile is $\mu - 1.28\,\sigma$, so a residual pedestal of $1.28\,\sigma$ per pixel is integrated over the whole frame; with no integration area, it pulls the centroid towards the frame center and inflates the widths. Simulated 1024×1280 frames (offset 64 ADU, read noise 5 ADU, shot noise):

| Spot | Quantity | True | Rainer (p10) | ISO 11146 |
|---|---|---|---|---|
| Far field, Gaussian $\sigma$ = 25 pix at (400, 300) | centroid | (400, 300) | (525.6, 410.7) | (400.00, 300.02) |
| | $\sigma_x$, $\sigma_y$ / $d_\sigma$ | 25, 25 / 100 | 293, 239 | $d_\sigma$ = 99.86 |
| Near field, top-hat R = 150 pix at (700, 550) | centroid | (700, 550) | (695.5, 547.2) | (699.99, 550.00) |
| | $\sigma_x$, $\sigma_y$ / $d_\sigma$ | 75, 75 / 300 | 125, 109 | $d_\sigma$ = 299.79 |

The bias grows with the ratio between frame area and spot flux: the method is adequate as a real-time pointing aid, not for quantitative measurements.

#### Gaussian fit (`lib/fitscmp_fit.py`)

Model: elliptical 2D Gaussian plus a constant background,

$$I(x, y) = b + A \exp\!\left[-\tfrac{1}{2}\left(\frac{u^2}{\sigma_1^2} + \frac{v^2}{\sigma_2^2}\right)\right], \quad u, v = \text{coordinates rotated by } \theta \text{ around } (x_0, y_0),$$

fitted by least squares (`scipy.optimize.least_squares`, trust-region reflective with bounds $\sigma > 0.3$ pix, $A \ge 0$) on the bounding box of the ISO integration area, starting from the ISO results. NaN pixels and **saturated pixels** (value $\ge$ `saturation_level`) are excluded, so the fit recovers the profile of a saturated spot from its unsaturated wings. Reported quantities: **fit center** $(x_0, y_0)$ (the center of the model, not the centroid), $d_{\sigma x}, d_{\sigma y} = 4\sigma$ and FWHM $= 2\sqrt{2\ln 2}\,\sigma$ along the principal axes, in the ISO convention for x, y and $\varphi$; $d_\sigma$, ellipticity and round flag as for ISO; fit amplitude $A$ and background $b$; total flux $2\pi A \sigma_1 \sigma_2$; number of fitted pixels.

**Uncertainties** (1σ) on the center and the widths come from the heteroscedasticity-consistent ("sandwich") covariance $(J^T J)^{-1} J^T \mathrm{diag}(r^2) J (J^T J)^{-1}$ [White 1980], with Jacobian $J$ and residuals $r$. The usual covariance, scaled by the mean residual variance, underestimates them by a factor 1.25–2.9 (larger at high signal), because shot noise makes the noise much larger on the spot than on the background; the sandwich form needs no noise model (the camera gain is not required). With it, the pulls (error / quoted uncertainty) have standard deviation 0.82–1.15 over 60 realizations at peak / read noise = 10, 40 and 400.

**Fit residual / noise** is the rms of the residuals divided by the background rms of the ISO analysis. It is 1.5–2 for bright Gaussian spots (shot noise on the spot), and of order tens for non-Gaussian profiles: 49 on the `NEARFIELD` top-hat, where the fit is not meaningful (fitted $d_\sigma$ = 231 instead of 180 pix). Use ISO 11146 for near fields.

#### Validation and accuracy

`tests/test_spot.py` (15 tests, about 40 s) checks, on simulated 512×640 frames: round and elliptical Gaussians (widths within 0.5–1%, azimuth within 0.5°, both orientation conventions of $\varphi$), a top-hat ($d_\sigma = 2R$ within 1%), sigma-clipped background with hot pixels outside the area, the in-frame flag, saturated plateaus, NaN pixels, the error on a constant frame, the exact equivalence of the Rainer method with the Acquisition Manager code, and for the Gaussian fit: widths, azimuth and flux of elliptical spots, exclusion of saturated pixels, robustness to hot pixels, calibration of the uncertainties (pulls), and the residual ratio of a non-Gaussian profile. On `test_spot.fits`, the spot moved by (+4.5, −2.8) pix and elongated to $\sigma$ = 24 × 18 pix at 30° gives $\Delta$ centroid = (+4.499, −2.813) pix, $d_{\sigma x}$ = 96.1, $d_{\sigma y}$ = 72.1 pix, $\varphi$ = 29.7°, flux ratio 1.080 (expected 1.08).

Second-moment widths are sensitive to noise in the integration area. Gaussian spot with $\sigma$ = 25 pix, 30 realizations:

| Peak / read noise | ISO $d_\sigma$ (mean ± std) | ISO centroid std | Fit $d_\sigma$ (mean ± std) |
|---|---|---|---|
| 10 | +2.8% ± 7.5% | 0.62 pix | −0.04% ± 0.24% |
| 40 | +0.1% ± 1.5% | 0.16 pix | −0.02% ± 0.09% |
| 400 | +0.02% ± 0.16% | 0.017 pix | 0.00% ± 0.02% |

For Gaussian spots the fit is 8 to 30 times less noisy (the gain is largest at low signal), because it uses the whole profile with the model as a constraint, while the second moments weight the noisy wings with $r^2$.

They are also sensitive to outliers far from the centroid: 30 hot pixels at random positions (1e-4 of a 512×640 frame, value 4095) bias $d_\sigma$ by +2.5% (median, max +5.3%) when they fall inside the integration area, while the Gaussian fit is unaffected (< 0.1%). ISO 11146 assumes a dark-corrected frame. Saturation flattens the profile and inflates the second-moment widths (80 → 82 pix on the `SATURATED` test spot, 100 → 107 pix in `tests/test_spot.py`); the fit, which excludes saturated pixels, gives 80.0 and 100.0 pix and recovers the amplitude within 0.1%.

**Integration area and frame size.** For a top-hat near field the 3 × widths rule is conservative: on `NEARFIELD` ($d_\sigma$ = 180 pix) the area is 540 pix wide and exceeds the 480-pix frame, so *Integration area in frame* is *no*, although the result is accurate (180.04 pix) because the background is flat. The flag reports non-conformity with the standard, not necessarily a wrong value.

Computation time: 0.7–1.3 s per 1024×1280 frame for ISO 11146; the Gaussian fit adds 0.3–0.8 s (it runs on the integration area only).

### 8.6 Frame calibration (`lib/fitscmp_calib.py`)

Applied to every image right after loading, before any other computation, in both views (*Calibration* frame of the GUI). Order: saturation mask on the raw values → master dark subtraction → hot-pixel correction.

#### Master dark

One master dark (for example the mean or median of several dark frames acquired with the same exposure time, gain and temperature as the images) is subtracted from all images: $I_\text{cal} = I - D$. It is read from the first 2D image of the selected FITS file (`dark_frame` in the configuration, or *Master dark...* in the GUI), must have the same shape as the images, and is not rescaled. After subtraction the background is close to zero; the spot analysis still estimates and removes any residual level.

The dark removes what a constant background cannot model: warm pixels and the **fixed pattern** of the detector. A fixed pattern well below the read noise already biases the second-moment widths, because the background model of ISO 11146 is a constant. Simulated Gaussian spot ($\sigma$ = 20 pix, read noise 5 ADU) with a sinusoidal pattern of period 251 pix, 10 realizations:

| Pattern amplitude | ISO $d_\sigma$ bias |
|---|---|
| 0 ADU | −0.01% ± 0.04% |
| 1 ADU | −1.78% ± 0.04% |
| 3 ADU | −4.42% ± 0.03% |

On `FARFIELD_HOT` (same pattern with 3 ADU amplitude and 40 warm pixels, [§12](#12-synthetic-test-data-make_test_fitspy)) ISO 11146 gives $d_\sigma$ = 76.9 pix on the raw frame, 76.4 pix with the hot-pixel filter only, and 79.93 pix after dark subtraction (true value 80). The Gaussian fit is much less sensitive (80.13 → 80.00 pix).

#### Saturation on raw values

The saturation mask is computed on the raw values, **before** the dark subtraction, and passed to the spot analysis: after subtraction a pixel at 4095 ADU would read about 4031 ADU and no longer be recognized. On the `SATURATED` test spot with the master dark subtracted, the 996 saturated pixels are still counted (and excluded from the Gaussian fit).

#### Hot-pixel filter

Optional correction of isolated hot pixels (*Hot-pixel filter*), useful when no dark is available or for residual defects. A pixel is corrected if:

1. it exceeds the **maximum** of its 8 neighbors by more than $n\,\sigma_\text{loc}$ (`hot_pixel_nsigma`, default 5), and
2. all its neighbors stay below half of its excess over the local background (median of the 8 neighbors).

Condition 2 protects the peak of a sampled spot, whose neighbors are comparable to it. The local noise $\sigma_\text{loc}^2 = a + b \cdot \text{level}$ (read + shot noise) is estimated from the image itself: the image minus its neighbor median is binned in 20 quantiles of the local level, and a line is fitted to the squared MAD of each bin, so the camera gain is not needed. With a single global noise value, normal shot-noise fluctuations on a bright spot exceeded the threshold (5 false detections on a 2000 ADU spot). Corrected pixels are replaced by the median of their neighbors, their number is reported in the table (*Hot pixels corrected*), and they are removed from the saturation mask.

Validation (`tests/test_calib.py` and simulations with read noise 5 ADU):

| Test | Result |
|---|---|
| False detections, 20 frames each (327 kpix/frame): far field, sharp spot ($\sigma$ = 1.5 pix, peak 4000 ADU), saturated spot, top-hat | 0, 1 (a background noise pixel, not the spot), 0, 0 |
| Pure noise, 1024×1280 | 0 |
| Peak of sharp spots with $\sigma$ = 1, 1.5, 3 pix | never flagged |
| Detection on the background | 30% at 6×, 100% from 10× the read noise |
| Detection on the spot core (~2000 ADU) | 0% up to 20×, 20% at 40×, 100% at 80× the read noise |
| ISO $d_\sigma$ with 30 hot pixels, after the filter | bias from > +1% to < 0.2% |

On a bright spot the threshold grows with the shot noise: small hot pixels there cannot be separated from photon noise and are left untouched, which is the conservative choice.

---

## 9. Output figures

### 9.1 *Subtraction* tab (2×2)

| Panel | Content | Color scale |
|---|---|---|
| A | image A in original units | from `clip_percentiles`: common to A and B without normalization, independent otherwise |
| B | image B in original units, after the shift | as A |
| A − B | difference of the **normalized** images | symmetric, $\pm$ upper percentile of $\lvert D \rvert$ |
| Histogram | histogram of $D$ (200 bins, log counts) with the statistics of [§8.2](#82-real-space-difference) | — |

With a normalization other than `none`, the titles of the A − B panel and of the histogram report it, and the difference statistics are in normalized units. Panels A, B and A−B share zoom and pan.

**Pixel coordinates.** Images are shown with `origin="lower"`; x is the column and y the row index, 0-based, with integer values at pixel centers. FITS tools such as DS9 use 1-based coordinates: add 1 to compare.

### 9.2 *Fourier* tab (2×3)

| Panel | Content |
|---|---|
| log10 P(A), log10 P(B) | 2D power spectra, common color scale |
| log10 P(A)/P(B) | 2D log ratio, symmetric color scale |
| Radial power spectrum | $\bar{P}_A$, $\bar{P}_B$, $\bar{P}_{A-B}$ vs $f$ (log–log) |
| Radial ratio | $\bar{P}_A / \bar{P}_B$ with reference line at 1 |
| Fourier Ring Correlation | FRC vs $f$ (log $f$), with reference line at 0 |

The figure title reports the NaN mode, the taper width and the fraction of pixels excluded from the FFT (after inpainting). The three 2D spectra share zoom and pan; the image extent is set at the pixel edges, so each frequency pixel is centered on its value (the DC pixel on $f_x = f_y = 0$).

### 9.3 *Images* tab (single-image view)

1 to 4 images in original units, laid out as 1×1, 1×2, 1×3 or 2×2, each with its own color scale from `clip_percentiles`. Zoom and pan are shared when all images have the same shape.

### 9.4 Spot overlay

With *Overlay* checked, each image panel (single-image view, and A and B in the *Subtraction* tab) shows the spot measured with the selected method:

| Marker | ISO 11146 | Gaussian fit | Rainer (p10) |
|---|---|---|---|
| red **+** | centroid | fit center | centroid |
| white **×** | peak position | peak position | peak position |
| red ellipse | full axes $d_{\sigma x}$, $d_{\sigma y}$ at $\varphi$ (the $2\sigma$ ellipse) | same, from the fit | full axes FWHM x, FWHM y, along the frame axes |
| white dashed polygon | integration area | fit region (pixel edges) | — |

In the comparison, B is measured before the shift but displayed after it: its overlay is translated by the applied shift, so that it lies on the displayed spot. The overlay does not change the axis limits.

### 9.5 Exported files

**Export...** in the result and batch windows asks for a base file name and writes (`export_format` = `csv`, default, or `parquet`, Brotli-compressed, which needs `pandas` and `pyarrow`):

| File | Window | Content |
|---|---|---|
| `<base>_summary.csv` | result | side table in tidy form, one row per image and quantity: `section`, `quantity`, `unit`, `column` (A, B, `B − A`, or `#k`), `value` (number; empty if not numeric), `text` (flags such as yes/no, iteration counts, and `ratio B/A` for the flux ratio, whose value is then the ratio) |
| `<base>_profiles.csv` | comparison | radial profiles: `f_cycles_per_pix`, `power_A`, `power_B`, `power_A_minus_B`, `ratio_A_over_B`, `frc` |
| `<base>_meta.json` | result | provenance and full-precision results: mode, image labels, all settings (`Settings.describe()`, the dark identified by its label), applied shift, pixel scale, hot-pixel counts, difference statistics, FFT excluded fraction and FRC 0.5 crossing, complete spot metrics of each image, export time and file list |
| `<base>_batch.csv` | batch | the batch table, with machine-friendly column names including the unit (e.g. `centroid_x_pix` or `centroid_x_um`, `frc_f50_cycles_per_pix`) |
| `<base>_batch_meta.json` | batch | reference, settings, number of processes and all batch rows at full precision |

Values are those shown in the table (in µm with *Units µm*), at full precision; NaN are written as empty cells in CSV and `null` in JSON. The figures can be saved as images from the Matplotlib toolbar.

---

## 10. Interpretation notes

- **Power spectra are shift-invariant**, while subtraction and FRC depend on phases: a pure translation leaves $P_A/P_B \simeq 1$ but lowers Pearson and FRC. Comparing both diagnostics separates *amplitude* differences (e.g. blur, smoothing, noise level) from *registration* differences.
- **After registration**, residual differences in subtraction and FRC reflect content rather than position. The automatic estimate aligns the cross-correlation peak, so a genuine asymmetric difference between the images (e.g. an asymmetric PSF or a moving feature) is partly absorbed into the shift: compare with the unregistered result, or with a manual shift, when this matters.
- **Additive white noise** in B produces a flat $\bar{P}_{A-B}$ and an FRC that drops where the signal power falls below the noise level.
- **Loss of resolution** (e.g. blur in B) shows up as $\bar{P}_A/\bar{P}_B > 1$ increasing with frequency.
- The lowest-frequency annuli contain few pixels and have a large variance; with the Hann window they are also slightly biased low by the mean subtraction.

---

## 11. API reference

```python
from lib.fitscmp_io import list_images, load_image, read_headers
from lib.fitscmp_analysis import difference, fourier_compare
from lib.fitscmp_registration import apply_shift, estimate_shift

refs = list_images("test_mef.fits")                 # -> list[ImageRef]
a, b = load_image(refs[2]), load_image(refs[5])     # -> float64 arrays (CUBE plane 0, SUBSHIFT)
dy, dx = estimate_shift(a, b)                       # shift to apply to B
b = apply_shift(b, dy, dx)                          # NaN on border strips
a, b, d, stats = difference(a, b, mode="median")    # normalized A, B, A-B, dict
fc = fourier_compare(a, b, window=True, nbins=64)   # dict (see below)
```

### `lib.fitscmp_io`

| Function | Returns |
|---|---|
| `list_images(path, expand_cubes=True)` | `list[ImageRef]` |
| `load_image(ref)` | `np.ndarray`, `float64`, shape `ref.shape` |
| `read_headers(path)` | `list[tuple[str, str]]` |

### `lib.fitscmp_analysis`

| Function | Returns |
|---|---|
| `image_stats(img)` | dict of min, max, mean, median, std, sum, `n_valid`, `n_nan` on the finite pixels (only `n_valid`, `n_nan` for an all-NaN image) |
| `normalize(img, mode="none", mask=None)` | normalized image, statistic on `img[mask]` if a mask is given; `ValueError` on unknown mode or invalid factor |
| `difference(a, b, mode="none")` | `(a_norm, b_norm, a_norm - b_norm, stats)`; `ValueError` on shape mismatch or no common finite pixel |
| `fourier_compare(a, b, window=True, nbins=64, nan_mode="interp", taper=16.0, interp_sigma=1.0)` | dict; `ValueError` on shape mismatch, unknown NaN mode or all-zero weights |

### `lib.fitscmp_mask`

| Function | Returns |
|---|---|
| `fft_weights(mask, window=True, taper=0.0)` | weight map $W$ (`float64`, shape of `mask`) |
| `inpaint(img, sigma=1.0)` | copy with NaN pixels filled; holes larger than the kernel stay NaN |
| `inpaint_common(a, b, sigma=1.0)` | `(a_filled, b_filled)`; pixels invalid in either image are filled in both, large holes stay NaN |

`NAN_MODES` is the tuple of accepted NaN modes.

### `lib.fitscmp_spot`

| Function | Returns |
|---|---|
| `spot_metrics(img, method="iso", saturation=None)` | dict of metrics (dispatches to the two functions below) |
| `iso_spot(img, saturation=None, max_iter=50)` | `background`, `background_rms`, `centroid_x/y`, `total_flux`, `d_sigma_x/y`, `d_sigma`, `phi_deg`, `ellipticity`, `round`, `area` (sides and angle), `area_in_frame`, `n_iter`, `converged`, `n_nan_in_area`, `peak_value`, `peak_x/y`, `n_peak`, `n_saturated` (if a level is given) |
| `rainer_spot(img, saturation=None)` | `background`, `centroid_x/y`, `total_flux`, `fwhm_x/y`, peak and saturation keys as above |
| `saturated_mask(img, saturation)` | boolean mask from a level, or the mask itself (`saturation` may be `None`, a level or a boolean mask in all spot functions) |
| `iso_widths(sxx, syy, sxy)` | `(d_sigma_x, d_sigma_y, phi)` from the second moments, $\varphi$ in radians |

Both raise `SpotError` (a `ValueError`) when the spot cannot be measured. `SPOT_METHODS` maps method ids to GUI labels; `ROUND_LIMIT` = 0.87.

### `lib.fitscmp_fit`

| Function | Returns |
|---|---|
| `gauss_spot(img, saturation=None)` | `background`, `amplitude`, `center_x/y` and `center_x/y_err`, `d_sigma_x/y` and `d_sigma_x/y_err`, `d_sigma`, `fwhm_x/y`, `phi_deg`, `ellipticity`, `round`, `total_flux`, `residual_rms`, `residual_ratio`, `n_fit`, `converged`, `fit_box` (x0, x1, y0, y1), peak and saturation keys |

### `lib.fitscmp_tables`

| Function | Returns |
|---|---|
| `measure_spots(imgs, method="iso", saturation=None)` | `[(metrics or None, error or None)]` per image; `saturation` is a level, or a list with one level or mask per image; metrics are computed once and shared by the table and the overlay |
| `batch_columns(method, saturation=None, hot_filter=False, scale=None)`, `batch_values(rows, cols, scale=None)` | batch columns (key, label, export name, is_length) and table values |
| `build_table(imgs, names, spots, method="iso", saturation=None, delta=False, scale=None, extra_stats=())` | `(rows, note)` for the side panel; rows are `(SECTION, title)` or `(label, values)`, with a B − A value when `delta=True`; lengths and positions in µm if `scale` [µm/pix] is given; `extra_stats` adds rows to the statistics section |

### `lib.fitscmp_pipeline`

```python
from lib.fitscmp_io import list_images
from lib.fitscmp_pipeline import Settings, run_batch

if __name__ == "__main__":  # required by the process pool (spawn)
    refs = list_images("test_spot.fits")
    rows = run_batch(refs[0], refs[1:], Settings(spot_method="gauss", auto_shift=True), workers=4)
```

| Function / class | Returns |
|---|---|
| `Settings(...)` | dataclass of all processing options (normalization, Fourier, spot method, saturation, dark, hot-pixel filter, registration); `describe()` gives a JSON-friendly dict |
| `load(refs, s)` | `(images, saturation masks, hot-pixel counts)`: loaded and calibrated |
| `run_single(refs, s)` | dict: `imgs`, `sats`, `n_hot`, `spots` |
| `run_comparison(pair, s)` | dict: `a`, `b` (calibrated, as recorded), `b_shifted`, `dy`, `dx`, `d`, `stats`, `fc`, `sats`, `n_hot`, `spots` |
| `compare_images(a, b, s)` | registration, difference and Fourier comparison of two calibrated images |
| `run_batch(reference, refs, s, workers=1, progress=None, cancel=None)` | list of row dicts sorted by index (0 = reference); `progress(done, total)`, `cancel()` → `True` raises `Cancelled` |
| `batch_row(k, ref, a, a_spot, s)` | one batch row |
| `frc_crossing(f, frc, level=0.5)` | first frequency where the FRC drops below `level` |

### `lib.fitscmp_export`

| Function | Returns |
|---|---|
| `export_result(base, table, meta, fc=None, fmt="csv")` | written paths: summary, profiles (if `fc`), metadata |
| `export_batch(base, columns, rows, meta, fmt="csv")` | written paths: batch table, metadata |
| `write_table(base, columns, rows, fmt="csv")`, `write_json(path, obj)`, `summary_rows(table)`, `to_jsonable(obj)` | helpers |

### `lib.fitscmp_calib`

| Function | Returns |
|---|---|
| `load_dark(path)` | `(array, label)` of the first 2D image of the file |
| `calibrate(img, dark=None, hot_filter=False, nsigma=5.0, saturation=None)` | `(calibrated image, saturation mask on raw values or None, number of corrected hot pixels)`; `ValueError` if the dark shape differs |
| `hot_pixels(img, nsigma=5.0)` | boolean mask of isolated hot pixels |

### `lib.fitscmp_windows`

`PairDialog`, `SelectDialog`, `BatchDialog`, `HeaderWindow`, `BatchWindow(parent, header, columns, values, fig, note="", export=None)` and `ResultWindow(parent, title, header, figs, table=None, readout=None, scale=None, export=None)`; `export` is a callable `export(base) → paths` run by the *Export...* button; `readout` is a list of `(name, array)` read under the cursor, and `ResultWindow.readout_text(event)` returns the readout for a Matplotlib mouse event (`None` outside image panels).

### `lib.fitscmp_plots`

| Function | Returns |
|---|---|
| `spot_overlay(m, dx=0.0, dy=0.0)` | overlay geometry (center, peak, ellipse, area) from spot metrics, translated by (dx, dy); `None` if `m` is `None` |
| `plot_subtraction(a, b, d, stats, cfg, norm="none", overlays=(None, None))` | `Figure` of [§9.1](#91-subtraction-tab-22) (A, B in original units, D normalized) |
| `plot_fourier(fc, cfg, note="")` | `Figure` of [§9.2](#92-fourier-tab-23) |
| `plot_single(imgs, titles, cfg, overlays=None)` | `Figure` of [§9.3](#93-images-tab-single-image-view), 1 to 4 images |
| `plot_batch(rows, method, scale=None)` | batch trends ([§6.6](#66-batch-comparison)) |

`NORMALIZATIONS` is the tuple of accepted normalization names.

Keys of the `fourier_compare` dictionary:

| Key | Shape | Content |
|---|---|---|
| `fx`, `fy` | `(nx,)`, `(ny,)` | centered frequency axes [cycles/pixel] |
| `logps_a`, `logps_b` | `(ny, nx)` | $\log_{10} P_A$, $\log_{10} P_B$ |
| `log_ratio` | `(ny, nx)` | $\log_{10}(P_A/P_B)$ |
| `f` | `(nbins,)` | radial bin centers |
| `prof_a`, `prof_b`, `prof_d` | `(nbins,)` | radial mean power of A, B, A−B |
| `frc` | `(nbins,)` | Fourier Ring Correlation |
| `excluded_frac` | scalar | fraction of pixels with no data in the FFT (after inpainting) |

### `lib.fitscmp_registration`

| Function | Returns |
|---|---|
| `estimate_shift(a, b, upsample=100, taper=16.0)` | `(dy, dx)` to apply to B; `ValueError` on shape mismatch or no valid pixel |
| `apply_shift(img, dy, dx, pad=16)` | shifted copy with NaN on border strips and around invalid pixels (the input itself for a zero shift); `ValueError` if the shift exceeds the image size |

Note: `difference()` expects **raw** images and returns normalized ones; `fourier_compare()` does not normalize, so it should receive the normalized images returned by `difference()` (as the GUI does).

---

## 12. Synthetic test data (`make_test_fits.py`)

Generates a noiseless 256×320 speckle-like pattern $s$ (random phase in a Gaussian pupil of width 0.05 cycles/pixel, intensity = $\lvert \text{field} \rvert^2$) scaled to a mean of $10^4$ (fixed seed, reproducible). Noise levels are fractions of $\sigma_s$. Files are written to the current directory.

| File | HDU | Content |
|---|---|---|
| `test_single.fits` | 0 | *a* = *s* + 1% noise, clipped at $6\times10^4$ (saturation) and stored as `uint16` (on disk: 16-bit integer with `BZERO = 32768`) |
| `test_mef.fits` | 0 | empty primary |
| | 1 `NOISY` | *b* = *a* + 5% white noise |
| | 2 `SHIFTED` | *a* circularly rolled by (+3, −2) pixels |
| | 3 `CUBE` | 2-plane cube [*a*, *b*] |
| | 4 `MASKED` | *b* with NaN defects: disk of radius 20 pix centered at (x, y) = (200, 128), column x = 60, 0.5% random pixels |
| | 5 `SUBSHIFT` | *s* shifted by (Δy, Δx) = (+1.3, −2.7) pix with the Fourier shift theorem, + independent 1% noise |

All extensions are `float32`. Indicative results with `CUBE` plane 0 as A (`mean` normalization, default Fourier settings, no registration):

| B | Tests | Pearson | `rms_rel` | Expected behavior |
|---|---|---|---|---|
| `NOISY` (with A = `CUBE` plane 1) | null case | 1.000 | 0 | $D = 0$, FRC = 1 |
| `NOISY` | white noise | 0.9988 | 0.050 | `rms_rel` equals the injected 5%; flat $\bar{P}_{A-B}$; FRC ≈ 1 at low $f$, → $0.01/\sqrt{0.01^2 + 0.05^2} \approx 0.20$ where noise dominates (the 1% noise is shared by A and B) |
| `MASKED` | NaN handling | 0.9988 | 0.050 | `masked` = 2.3%, 0.8% of pixels excluded from the FFT; spectra and FRC close to the `NOISY` case |
| `SHIFTED` | integer shift | 0.52 | 0.98 | radial ratio ≈ 1, FRC strongly reduced |
| `SUBSHIFT` | sub-pixel shift | 0.65 | 0.84 | as above |
| `test_single` | integer scaling, saturation | 0.9985 | 0.055 | 153 saturated pixels, `max_abs` ≫ `rms` |

With *Auto* registration (upsample 100):

| B | True shift to apply | Estimate | Pearson | `rms_rel` | Notes |
|---|---|---|---|---|---|
| `NOISY`, `MASKED` | (0, 0) | (0.00, 0.00) | unchanged | unchanged | no spurious shift |
| `SHIFTED` | (−3, +2) | (−2.99, +2.00) | 1.0000 | 0.003 | shared noise realization; 1.8% of pixels in the border strips |
| `SUBSHIFT` | (−1.3, +2.7) | (−1.30, +2.70) | 0.9999 | 0.0155 | at the noise floor $\sqrt{2} \times 1\%$ ≈ 0.014 |
| `SUBSHIFT` + `MASKED` defects | (−1.3, +2.7) | (−1.30, +2.68) | 0.9998 | 0.018 | 5.6% of pixels masked (defects, growth, strips) |

`test_spot.fits` contains laser spots on a 480×640 camera frame (dark offset 64 ADU, read noise 5 ADU, shot noise with gain 2 e⁻/ADU, rounded and clipped to 12 bits, `uint16`):

| HDU | Content | ISO 11146 result |
|---|---|---|
| 1 `FARFIELD` | round Gaussian, $\sigma$ = 20 pix, peak 2500 ADU, at (300, 240) | centroid (299.99, 240.02), $d_\sigma$ = 80.0, round |
| 2 `FARFIELD_MOVED` | Gaussian $\sigma$ = 24 × 18 pix, major axis at 30°, at (304.5, 237.2) | centroid (304.48, 237.20), $d_{\sigma x}$ = 96.1, $d_{\sigma y}$ = 72.1, $\varphi$ = 29.7°, $\varepsilon$ = 0.750 |
| 3 `NEARFIELD` | uniform disk of radius 90 pix at (320, 240), 1800 ADU | centroid (319.99, 240.00), $d_\sigma$ = 180.1 (= 2R) |
| 4 `SATURATED` | as `FARFIELD` with peak 6000 ADU, clipped at 4095 | 996 saturated pixels, $d_\sigma$ = 82.2 (inflated by the plateau) |
| 5 `FARFIELD_HOT` | as `FARFIELD`, on the dark pattern of `test_dark.fits` (offset 64 ADU + 3 ADU sinusoidal fixed pattern + 40 warm pixels of 200–1500 ADU) | $d_\sigma$ = 76.9 raw, 76.4 with the hot-pixel filter only, 79.93 after dark subtraction |

`test_dark.fits` is the master dark of that pattern (mean of 16 frames with read noise, `float32`, `NCOMBINE` = 16).

---

## 13. Known limitations

**Method**
- Registration handles translations only (no rotation, scale or distortion); the estimate is quantized to $1/u$ pixel.
- Fourier resampling assumes band-limited data: sharp features and the edges of large filled holes produce ringing; taking the real part after the shift slightly alters the Nyquist component of even-sized axes.
- Images must have identical shape; no cropping, binning or resampling. The pair dialog does not filter incompatible images.
- Masked spectra are pseudo-spectra (no mode-coupling deconvolution); residual leakage from large holes remains at high frequency, especially for spectra with a large dynamic range. The apodization also reduces the effective area around each hole.
- `interp` mode replaces invalid pixels with interpolated values in the Fourier analysis: appropriate for isolated defects, not for extended masked regions (these are larger than the kernel and remain masked).
- Frequencies $f \ge 0.5$ cycles/pixel (corners and the exact Nyquist ring) are excluded from radial profiles.
- No uncertainties on radial profiles and no FRC threshold curves (e.g. 1/7 or half-bit criteria).
- When A−B is identically zero, its radial profile is not visible on the log–log plot.

**Software**
- Each image is held as `float64` plus several complex arrays of the same size; each batch process holds a copy of the reference image.
- The process pool uses the *spawn* start method: when the pipeline is called from a script, `run_batch(..., workers > 1)` must be called under `if __name__ == "__main__":` and the script must be a file (not code read from standard input), as usual with `multiprocessing`.
- A single comparison cannot be interrupted: *Cancel* discards its result when it ends. Only one job runs at a time.
- Figures are not exported by *Export...* (use the Matplotlib toolbar).
- Configuration values other than `normalization` and `nan_mode` are not validated. A non-numeric or negative bin count typed in the GUI produces an error dialog; 0 produces empty radial plots.
- Arrays with more than three dimensions are ignored.
- The Gaussian fit models a single elliptical Gaussian: it is not meaningful for top-hat near fields or multimode speckled profiles (check *Fit residual / noise*).
- *Units µm* assumes square pixels and applies to the table only.
- The master dark is subtracted as is (no scaling with exposure time), and must match the images in shape, exposure, gain and temperature.
- The hot-pixel filter corrects isolated pixels only (adjacent pairs or clusters are not detected) and leaves small hot pixels on bright spots untouched ([§8.6](#86-frame-calibration-libfitscmp_calibpy)).
- Spot analysis assumes a single spot per image and a dark-corrected frame; hot pixels inside the integration area and saturation bias the second-moment widths, and low signal-to-noise ratios increase their scatter ([§8.5](#85-spot-analysis-libfitscmp_spotpy)).
- The single-image view shows at most 4 images; zoom and pan are not shared between images of different shapes.

---

## 14. References

- White, H. (1980). A heteroskedasticity-consistent covariance matrix estimator and a direct test for heteroskedasticity. *Econometrica*, 48(4), 817–838. doi:10.2307/1912934
- Harris, F. J. (1978). On the use of windows for harmonic analysis with the discrete Fourier transform. *Proc. IEEE*, 66(1), 51–83. doi:10.1109/PROC.1978.10837
- Saxton, W. O., & Baumeister, W. (1982). The correlation averaging of a regularly arranged bacterial cell envelope protein. *J. Microscopy*, 127(2), 127–138. doi:10.1111/j.1365-2818.1982.tb00405.x
- van Heel, M., & Schatz, M. (2005). Fourier shell correlation threshold criteria. *J. Struct. Biol.*, 151(3), 250–262. doi:10.1016/j.jsb.2005.05.009
- ISO 11146-1:2021. *Lasers and laser-related equipment — Test methods for laser beam widths, divergence angles and beam propagation ratios — Part 1: Stigmatic and simple astigmatic beams.* International Organization for Standardization.
- ISO 11146-3:2004. *Lasers and laser-related equipment — Test methods for laser beam widths, divergence angles and beam propagation ratios — Part 3: Intrinsic and geometrical laser beam classification, propagation and details of test methods.* International Organization for Standardization.
- Guizar-Sicairos, M., Thurman, S. T., & Fienup, J. R. (2008). Efficient subpixel image registration algorithms. *Optics Letters*, 33(2), 156–158. doi:10.1364/OL.33.000156
- Hivon, E., Górski, K. M., Netterfield, C. B., et al. (2002). MASTER of the cosmic microwave background anisotropy power spectrum. *ApJ*, 567(1), 2–17. doi:10.1086/338126
- Astropy documentation — Convolution and NaN interpolation (`interpolate_replace_nans`): <https://docs.astropy.org/en/stable/convolution/>
- Astropy documentation — *Image data* and scaled data handling: <https://docs.astropy.org/en/stable/io/fits/usage/image.html>
- scikit-image documentation — `skimage.registration.phase_cross_correlation`: <https://scikit-image.org/docs/stable/api/skimage.registration.html>
- NumPy documentation — Discrete Fourier Transform (`numpy.fft`): <https://numpy.org/doc/stable/reference/routines.fft.html>
