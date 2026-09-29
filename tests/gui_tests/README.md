# fitscmp GUI tests (not part of the repository)

Ad-hoc scripts used during development to validate the Tk GUI. They drive `fitscmp_gui.App`
programmatically (dialogs closed with `after()` callbacks, private attributes such as `app._job`,
`dlg._cb`), so they are fragile by design and meant to be ported to a proper `tests/gui/` suite.

Run from the project root, after `python make_test_fits.py`, with a virtual X server:

    xvfb-run -a python /path/to/gui_testN.py

| Script | Checks |
|---|---|
| gui_test.py | swap A/B with two images, dialog defaults with swap |
| gui_test2.py | NaN modes (interp/mask), Fourier figure title |
| gui_test3.py | registration: auto, manual edit, reset |
| gui_test4.py | single-image view, raw values, max 4 images, raw A/B panels |
| gui_test5.py | spot table: ISO deltas, Rainer rows, saturation row |
| gui_test6.py | Gaussian fit rows, overlay (also shifted B), µm conversion, overlay off |
| gui_test7.py | master dark, saturation on raw values, hot-pixel filter, dark shape error, pixel readout and pins |
| gui_test8.py | batch (1 and 2 processes, identical), drift recovery, CSV/Parquet export, busy state, cancel (needs `if __name__ == "__main__"`) |

Requirements beyond requirements.txt: python3-tk, xvfb (xvfb-run), pandas + pyarrow for the Parquet part of gui_test8.
