# Changelog

## [Unreleased]

### Step 8 — export, background execution, batch comparison

#### Added
- `lib/fitscmp_pipeline.py`: GUI-independent processing (`Settings`, `run_single`, `run_comparison`, `run_batch`, `frc_crossing`); the GUI only collects the options and displays the results.
- Background execution: jobs run in a worker thread (Tk used only from the GUI thread), progress bar and *Cancel* in the status line, one job at a time.
- Batch comparison (*Batch...*): every image compared with a reference; table with position, displacement, widths, flux ratio, registration shift, rms_rel, Pearson and FRC 0.5 crossing, and trend figure; process pool (`batch_workers`, spawn start method) with results identical to the sequential run.
- `lib/fitscmp_export.py` and *Export...* buttons: tidy summary table, radial profiles and batch table as CSV or Brotli-compressed Parquet (`export_format`), JSON metadata with settings and full-precision results.
- `tests/test_pipeline.py` (5 tests: batch drift recovery, pool vs sequential, cancel, FRC crossing, export round trip).

#### Changed
- *Auto* registration label: "estimate per comparison" (the estimate is also made for each batch image).
- NaN values are shown as a dash in all tables.

### Step 7 — calibration, pixel readout, windows module

#### Added
- `lib/fitscmp_calib.py`: master dark subtraction (`dark_frame` key, *Master dark...* in the GUI) and hot-pixel filter (isolated pixels above the maximum of their neighbors by `hot_pixel_nsigma` times a read + shot noise model estimated from the image); *Calibration* frame applied to all images.
- Saturation evaluated on the raw values: spot functions accept a saturation level or a boolean mask (`saturated_mask()`); per-image masks in `measure_spots()`.
- Pixel readout above the result tabs: coordinates and values of all images at the pixel under the cursor, left click to pin (up to 5).
- *Hot pixels corrected* row in the table; calibration reported in the window header.
- `make_test_fits.py`: `test_dark.fits` (master dark) and the `FARFIELD_HOT` spot (generated last: previous test data unchanged).
- `tests/test_calib.py` (6 tests).

#### Changed
- Dialogs, header browser and result window moved to `lib/fitscmp_windows.py` (`fitscmp_gui.py`: 434 → 305 lines, including the new calibration controls).
- GUI frames ordered as the processing: Options, Calibration, Masked pixels (FFT), Spot analysis, Registration.

### Step 6 — Gaussian fit, overlay, µm units

#### Added
- `lib/fitscmp_fit.py`: elliptical 2D Gaussian + background fit on the ISO integration area, saturated pixels excluded, heteroscedasticity-consistent (sandwich) uncertainties, residual / noise ratio; *Gaussian fit* option of the spot method.
- Spot overlay (center, peak position, ellipse, integration area or fit region) on the image panels, *Overlay* option; B's overlay follows the registration shift.
- *Units µm* option and `pixel_scale` configuration key: positions and widths in the table in µm.
- Fit tests in `tests/test_spot.py` (15 tests in total).

#### Changed
- `lib/fitscmp_tables.py`: spots are measured once (`measure_spots()`) and shared by table and overlay; `build_table()` takes the measured spots and an optional pixel scale.

### Step 5 — spot metrics (ISO 11146, Rainer p10)

#### Added
- `lib/fitscmp_spot.py`: ISO 11146 spot analysis (first and second moments, widths dσx, dσy, dσ along the principal axes, azimuth, ellipticity and round flag, iterative background and integration area of 3 x widths, in-frame flag) and the Rainer (p10) method of the Acquisition Manager (bit-identical); peak value and position, total flux, saturated pixels.
- `lib/fitscmp_tables.py`: side-panel table with *Image statistics* and *Spot* sections, B − A column in the comparison (ratio for the flux), definitions and errors in a note.
- GUI frame *Spot analysis* (method); YAML keys `spot_method`, `saturation_level`.
- `tests/test_spot.py`: validation on simulated spots.
- `make_test_fits.py`: `test_spot.fits` with far-field, moved/elliptical, near-field and saturated spots on a 12-bit camera frame.

#### Changed
- The comparison table refers to the images as recorded (B before the registration shift).

### Step 4 — single-image view, original units

#### Added
- *Single images* option: 1 to 4 images shown side by side in original units (`SelectDialog`), used automatically when only one image is available; the run button reads *Show* in this mode.
- Side table in the result window with the statistics of each image in original units (`image_stats()`).

#### Changed
- Comparison: panels A and B show the images in original units (B after the shift); only A − B and its statistics use the normalization, which is reported in the titles.
- Plotting functions moved from `fitscmp_gui.py` to `lib/fitscmp_plots.py`; new `plot_single()`.
- `ResultWindow` takes a window title, a header text and an optional table.

#### Fixed
- On Linux screens with a DPI different from 96, figures could be drawn larger than their canvas and clipped on the right (Matplotlib TkAgg device-pixel-ratio ordering); each figure is now resized to its widget when mapped.

### Step 3 — registration

#### Added
- `lib/fitscmp_registration.py`: `estimate_shift()` (cross-correlation with upsampled-DFT refinement, on Hann-windowed images with the common apodized mask) and `apply_shift()` (Fourier shift theorem with mirror padding; border strips and pixels near NaN input set to NaN).
- GUI frame *Registration*: dy, dx fields, *Auto* (estimate at Compare and fill the fields), *Reset*; the applied shift is shown in the result window and in the status bar. YAML key `registration_upsample`.
- `lib/fitscmp_mask.inpaint()` for single images (`inpaint_common()` now uses it).
- scikit-image dependency.

#### Fixed
- `docs/DOCUMENTATION.md`: a leftover fragment of the old test-data section (introduced in step 2) has been removed.

### Step 2 — masked pixels in the Fourier analysis

#### Added
- `lib/fitscmp_mask.py`: common weight map for A, B and A−B (validity mask × cosine apodization × Hann window) and inpainting of pixels invalid in either image.
- `fourier_compare()`: arguments `nan_mode` (`interp` default, or `mask`), `taper` (default 16 pix) and `interp_sigma` (default 1 pix); output key `excluded_frac`.
- GUI frame *Masked pixels (FFT)* (mode, taper) and a figure title reporting the settings and the excluded fraction; new YAML keys `nan_mode`, `mask_taper`, `interp_sigma`.
- `masked` statistic: fraction of pixels not finite in both images.
- `make_test_fits.py`: extensions `MASKED` (NaN defects) and `SUBSHIFT` (sub-pixel shift, for the registration step), appended after the existing ones.
- SciPy dependency (distance transform).

#### Changed
- Power spectra are normalized by ΣW² with the full weight map (reduces to the previous normalization without NaN pixels).
- `make_test_fits.py`: the reference image now includes a 1% noise floor; without it the spectrum spans ~12 decades, which no masking method can preserve.

#### Removed
- Independent mean filling of NaN pixels in A and B before the FFT.

### Step 1 — statistics, GUI fixes, packaging

#### Changed
- `difference()`: normalization factors are computed on the common finite mask of A and B (`normalize()` gains an optional `mask` argument).
- `rms_rel` is now RMS(A−B) / std(A) instead of RMS(A−B) / RMS(A), so it no longer depends on the image pedestal.
- Statistics of constant images return NaN/inf without runtime warnings.

#### Added
- *Swap A/B* option in the GUI: inverts the order when exactly two images are found, and the default selection of the pair dialog otherwise.
- `requirements.txt`, `.gitignore`, `CHANGELOG.md`.
- Technical documentation (`docs/DOCUMENTATION.md`) and `README.md`.

#### Fixed
- 2D power spectra were displayed with a half-pixel offset (extent now set at the pixel edges).
