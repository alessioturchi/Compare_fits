"""FITS I/O helpers: 2D image discovery, loading and header extraction."""
from __future__ import annotations

import os
from dataclasses import dataclass

import numpy as np
from astropy.io import fits

_IMAGE_HDUS = (fits.PrimaryHDU, fits.ImageHDU, fits.CompImageHDU)


@dataclass(frozen=True)
class ImageRef:
    """Pointer to a 2D image inside a FITS file."""
    path: str
    hdu: int
    plane: int | None  # index along the first numpy axis for 3D cubes, None for 2D HDUs
    extname: str
    shape: tuple

    def label(self):
        s = f"{os.path.basename(self.path)} [HDU {self.hdu}"
        if self.extname:
            s += f" '{self.extname}'"
        if self.plane is not None:
            s += f", plane {self.plane}"
        return s + f"] {self.shape[0]}x{self.shape[1]}"


def list_images(path, expand_cubes=True):
    """List all 2D images in a FITS file (shape read from headers, no data loaded)."""
    refs = []
    with fits.open(path) as hdul:  # default memmap=None: lazy, scaling allowed
        for i, hdu in enumerate(hdul):
            if not isinstance(hdu, _IMAGE_HDUS):
                continue
            shape = tuple(hdu.shape)
            name = hdu.name if hdu.name not in ("", "PRIMARY") else ""
            if len(shape) == 2:
                refs.append(ImageRef(path, i, None, name, shape))
            elif len(shape) == 3 and expand_cubes:
                refs += [ImageRef(path, i, k, name, shape[1:]) for k in range(shape[0])]
    return refs


def load_image(ref):
    """Load a 2D image as native-endian float64 (BZERO/BSCALE applied, copied from the file)."""
    with fits.open(ref.path) as hdul:  # explicit memmap=True forbids BZERO scaling
        data = hdul[ref.hdu].data
        arr = data if ref.plane is None else data[ref.plane]
        return np.array(arr, dtype=np.float64)


def read_headers(path):
    """Return a list of (title, header text) for every HDU of a FITS file."""
    with fits.open(path) as hdul:
        return [(f"HDU {i}: {hdu.name} ({type(hdu).__name__})",
                 hdu.header.tostring(sep="\n", padding=False))
                for i, hdu in enumerate(hdul)]
