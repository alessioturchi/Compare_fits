"""Create synthetic FITS files to validate fitscmp_gui.py."""
import numpy as np
from astropy.io import fits

rng = np.random.default_rng(1)
ny, nx = 256, 320
fy, fx = np.meshgrid(np.fft.fftfreq(ny), np.fft.fftfreq(nx), indexing="ij")


def speckle(width=0.05):
    """Speckle-like intensity: random phase in a Gaussian pupil (width in cycles/pix), via FFT."""
    pupil = np.exp(-(fx**2 + fy**2) / (2 * width**2))
    field = np.fft.ifft2(pupil * np.exp(2j * np.pi * rng.random((ny, nx))))
    return np.abs(field) ** 2


def fourier_shift(img, dy, dx):
    """Circular sub-pixel shift by the Fourier shift theorem (amplitudes preserved)."""
    return np.real(np.fft.ifft2(np.fft.fft2(img) * np.exp(-2j * np.pi * (fy * dy + fx * dx))))


def with_defects(img):
    """Copy of img with NaN defects: a disk (r = 20 pix), a bad column and 0.5% bad pixels."""
    yy, xx = np.mgrid[:ny, :nx]
    out = img.copy()
    out[(yy - 128) ** 2 + (xx - 200) ** 2 < 20 ** 2] = np.nan
    out[:, 60] = np.nan
    out[rng.random((ny, nx)) < 0.005] = np.nan
    return out


s = speckle()
s = 1e4 * s / s.mean()                                   # noiseless pattern
noise = lambda f: rng.normal(0, f * s.std(), s.shape)    # white noise, fraction f of the pattern std
a = s + noise(0.01)                                      # reference: pattern + 1% noise floor
b = a + noise(0.05)                                      # same realization + 5% extra noise
c = np.roll(a, (3, -2), axis=(0, 1))                     # integer circular shift
d = with_defects(b)                                      # noisy pattern with NaN defects
e = fourier_shift(s, 1.3, -2.7) + noise(0.01)            # sub-pixel shift, independent 1% noise

# File 1: single 2D image stored as 16-bit integers with BZERO; clipping at 6e4 mimics saturation
fits.PrimaryHDU(np.clip(a, 0, 6e4).astype(np.uint16),
                header=fits.Header({"OBJECT": "A"})).writeto("test_single.fits", overwrite=True)
# File 2: multi-extension with 2D images and a 3D cube (new extensions appended at the end)
fits.HDUList([fits.PrimaryHDU(),
              fits.ImageHDU(b.astype(np.float32), name="NOISY"),
              fits.ImageHDU(c.astype(np.float32), name="SHIFTED"),
              fits.ImageHDU(np.stack([a, b]).astype(np.float32), name="CUBE"),
              fits.ImageHDU(d.astype(np.float32), name="MASKED"),
              fits.ImageHDU(e.astype(np.float32), name="SUBSHIFT")]).writeto("test_mef.fits", overwrite=True)
# File 3: laser spots on a 12-bit camera frame (dark offset, read and shot noise), uint16
sy, sx = 480, 640
yy, xx = np.mgrid[:sy, :sx].astype(float)


def camera(signal):
    """Dark offset 64 ADU, read noise 5 ADU, shot noise (gain 2 e-/ADU), 12-bit saturation."""
    frame = signal + 64 + rng.normal(0, 5, signal.shape) + rng.normal(0, np.sqrt(np.maximum(signal, 0) / 2))
    return np.clip(np.round(frame), 0, 4095).astype(np.uint16)


def gauss(x0, y0, s1, s2, theta, peak):
    """Elliptical Gaussian, major axis at theta [deg] from x."""
    t = np.radians(theta)
    u = (xx - x0) * np.cos(t) + (yy - y0) * np.sin(t)
    v = -(xx - x0) * np.sin(t) + (yy - y0) * np.cos(t)
    return peak * np.exp(-0.5 * ((u / s1) ** 2 + (v / s2) ** 2))


spots = fits.HDUList([fits.PrimaryHDU(),
                      fits.ImageHDU(camera(gauss(300.0, 240.0, 20, 20, 0, 2500)), name="FARFIELD"),
                      fits.ImageHDU(camera(gauss(304.5, 237.2, 24, 18, 30, 2500)), name="FARFIELD_MOVED"),
                      fits.ImageHDU(camera(1800.0 * (((xx - 320) ** 2 + (yy - 240) ** 2) < 90 ** 2)), name="NEARFIELD"),
                      fits.ImageHDU(camera(gauss(300.0, 240.0, 20, 20, 0, 6000)), name="SATURATED")])

# File 4: master dark (mean of 16 dark frames) with a fixed pattern and warm pixels; the
# FARFIELD_HOT spot contains the same pattern (generated last: the data above are unchanged)
pattern = 64.0 + 3.0 * np.sin(xx / 40.0)
warm_y, warm_x = rng.integers(0, sy, 40), rng.integers(0, sx, 40)
pattern[warm_y, warm_x] += rng.uniform(200.0, 1500.0, 40)
dark = np.mean([pattern + rng.normal(0, 5, pattern.shape) for _ in range(16)], axis=0)
fits.PrimaryHDU(dark.astype(np.float32), header=fits.Header({"OBJECT": "MASTER DARK", "NCOMBINE": 16})
                ).writeto("test_dark.fits", overwrite=True)
sig = gauss(300.0, 240.0, 20, 20, 0, 2500)
hot = sig + pattern + rng.normal(0, 5, sig.shape) + rng.normal(0, np.sqrt(sig / 2))
spots.append(fits.ImageHDU(np.clip(np.round(hot), 0, 4095).astype(np.uint16), name="FARFIELD_HOT"))
spots.writeto("test_spot.fits", overwrite=True)
print("Written test_single.fits (1 image), test_mef.fits (4 images + 2-plane cube), "
      "test_spot.fits (5 spots), test_dark.fits (master dark)")
