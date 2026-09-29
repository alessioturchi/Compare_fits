# fitscmp

A small Tkinter tool to view 2D images from FITS files and compare them, in real space and in Fourier space.

## Features

- Reads 2D images from any image HDU (`PrimaryHDU`, `ImageHDU`, tile-compressed `CompImageHDU`), including multi-extension files and individual planes of 3D cubes; `BSCALE`/`BZERO` are applied.
- Single-image view (1 to 4 images, automatic when only one image is available) with statistics in original units.
- Optional flux normalization (`mean`, `median`, `max`, `sum`), computed on the pixels valid in both images.
- **Real space**: A and B in original units, A−B map, histogram of the difference, RMS, relative RMS, maximum deviation, Pearson correlation.
- **Fourier space**: 2D log power spectra and their log ratio, radial power spectra of A, B and A−B, radial ratio, Fourier Ring Correlation (FRC). Optional Hann window.
- Consistent handling of NaN pixels in the Fourier analysis: common mask for A and B, inpainting of small defects, apodization of masked regions.
- Laser-spot metrics in a side table: ISO 11146 centroid, second-moment widths, azimuth and ellipticity, peak value and position, total flux, saturation; alternatively an elliptical Gaussian fit with uncertainties (saturated pixels excluded), or the Rainer (p10) method of the Acquisition Manager for comparison. Optional overlay on the images and values in µm.
- Frame calibration: master dark subtraction and an optional filter for isolated hot pixels; saturation is checked on the raw values.
- Pixel readout under the cursor with the values of all images at that pixel, and pinned readouts.
- Sub-pixel registration of B onto A: manual shift or automatic estimate (upsampled cross-correlation), applied with the Fourier shift theorem.
- Batch comparison of a series of images with a reference (table and trends of position, width, flux and difference), parallel processing.
- Export of the numerical results (CSV or Brotli-compressed Parquet) with JSON provenance metadata.
- Computations in the background: the interface stays responsive, with progress bar and cancel.
- FITS header browser.
- The I/O and analysis modules (`lib/`) are GUI-independent and can be used from scripts.

## Requirements

Python 3 with NumPy, SciPy, scikit-image, Astropy, Matplotlib (≥ 3.5), PyYAML and Tkinter.

```bash
pip install -r requirements.txt
```

On some Linux distributions Tkinter must be installed separately (e.g. `sudo apt install python3-tk`). Parquet export additionally needs `pandas` and `pyarrow` (optional).

## Usage

```bash
python fitscmp_gui.py
```

Run from the project root. Add one or more FITS files, press **Compare** and, if more than two images are available, choose images A and B (*Swap A/B* inverts the default order). Check *Single images* to view images without comparing them. Results open in a window with a *Subtraction* and a *Fourier* tab.

To try it on synthetic data:

```bash
python make_test_fits.py   # creates test_single.fits, test_mef.fits, test_spot.fits and test_dark.fits
python tests/test_spot.py  # validation of the spot metrics
python tests/test_calib.py # validation of dark subtraction and hot-pixel filter
python tests/test_pipeline.py  # validation of batch, process pool and export
```

Default options are set in `fitscmp.yaml`; normalization, windowing, number of radial bins, NaN handling and registration can also be changed in the GUI.

### Scripted use

```python
from lib.fitscmp_io import list_images, load_image
from lib.fitscmp_analysis import difference, fourier_compare
from lib.fitscmp_registration import apply_shift, estimate_shift

refs = list_images("test_mef.fits")
a, b = load_image(refs[2]), load_image(refs[5])
b = apply_shift(b, *estimate_shift(a, b))
a, b, d, stats = difference(a, b, mode="median")
fc = fourier_compare(a, b, window=True, nbins=64)
```

## Project layout

```
fitscmp_gui.py        GUI and plotting
fitscmp.yaml          default settings
make_test_fits.py     synthetic test data
lib/                  FITS I/O, analysis, spot metrics and plotting
tests/                validation tests
docs/                 documentation
requirements.txt      Python dependencies
CHANGELOG.md          history of changes
```

## Documentation

Methods, formulas, configuration, API and known limitations are described in [docs/DOCUMENTATION.md](docs/DOCUMENTATION.md).

## Limitations

The two images must have the same shape; registration handles translations only (no rotation or scale). Spot metrics assume one spot per image and a dark-corrected frame. See the documentation for the full list.

## License

<!-- TODO: choose a license and add the LICENSE file -->

## Author

Alessio Turchi
