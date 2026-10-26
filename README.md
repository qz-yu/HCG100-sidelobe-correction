# HCG100 sidelobe correction

Public validation code and input data for the FAST H I sidelobe correction in the HCG100 field. The workflow is based only on the three notebooks supplied with this release.

## Run

Use Python 3.12 or later in a separate environment. Numerical versions used for this release are also available in `requirements-validated.txt`:

```bash
python -m pip install -r requirements.txt
python prepare_data.py --download
python validate.py
jupyter lab
```

Open the notebooks from the repository root and run every cell in this order:

1. **01_prepare_vla_scaled_truth_maps.ipynb**: mask the VLA cube, convolve it to the beam in the FAST header, and interpolate/reproject it to the FAST velocity and spatial grid. Check the FAST/VLA peak ratios and construct scaled source maps at native VLA resolution.
2. **02_sidelobe_maps_5200_5600.ipynb**: convolve each native-resolution source map with each of the 19 real and Gaussian beam responses. Subtract the Gaussian response, reproject to the FAST grid, and assemble the map using the 32 pointing tables.
3. **03_sidelobe-subtraction-smooth-2fov-VLA.ipynb**: distribute each integrated sidelobe map across velocity using its smoothed VLA aperture spectrum, subtract the summed model from the FAST cube, and apply the original spatial smoothing kernel.

All paths are relative to the repository root. The notebooks write intermediate products into the existing data directories and into `sidelobe-result-fits/`; rerunning them overwrites those products. Input files are listed in `input_manifest.json`.

## Method and parameters

The reference maps use the masked VLA emission over **5200-5600 km/s**. For the beam convolution, the native VLA moment maps are converted from mJy per VLA beam km/s to mJy arcmin^-2 km/s using the VLA Gaussian beam area, `1.1331 * BMAJ * BMIN`.

| Source | RA (deg) | Dec (deg) | Aperture radius (arcmin) | Adopted scale factor |
| --- | ---: | ---: | ---: | ---: |
| HCG100 | 0.3022699 | 13.0891334 | 4.0 | 1.2846954129518 |
| NGC7810 | 0.5763757 | 12.9642823 | 2.2 | 1.1769393092739915 |
| MRK935 | 0.4504456 | 13.1028654 | 1.75 | 1.0799621913039583 |
| AGC105092 | 0.4095431 | 12.9043733 | 2.2 | 1.4756065799816787 |

The scale factors are the fixed values recorded in notebook 01; running its peak-ratio diagnostic does not automatically replace them. The diagnostic uses a circular pixel mask. The channel spectra in notebook 03 use fractional-pixel aperture photometry.

Both beam kernels are resampled to the VLA pixel scale and divided by the **real-beam peak**. Their array centering is used as stored; the implementation applies no CRPIX shift. The integrated model is

```text
S_source = pixel_area x [truth_source * real_beam - truth_source * Gaussian_beam]
```

where `*` denotes convolution and `x` denotes scalar multiplication. At each sampled FAST pixel, contributions from repeated pointings are averaged. Unsampled pixels receive zero before smoothing. The 19-beam sum saved as `diagnostic_sum19_*` is a diagnostic, not the assembled correction.

For each source, negative aperture-spectrum values are set to zero for normalization. The resulting weights satisfy `sum(weight * dv) = 1`, with `dv = 20 km/s` for the supplied cube. The source model is `S_source * weight(v)`; its spatial pattern is therefore fixed across channels. Nearest-channel matching requires a separation of less than half a FAST channel.

The final smoothing uses `Gaussian2DKernel(x_stddev=0.8344, y_stddev=0.9631)` in pixels, retains the original NaN mask, and labels the output beam as 4 arcmin x 4 arcmin. These are the settings in notebook 03.

## Included data

| Directory | Inputs |
| --- | --- |
| `VLA/` | HCG100 H I cube and its mask |
| `rawcube/` | FAST two-field cube before this correction |
| `moments/` | FAST moment-0 map used as the spatial header reference |
| `sidelobe-corr/beam-big/` | 19 real response arrays (`*-beam-kernal.fits`) and 19 Gaussian arrays (`*-beam-clean2.fits`) |
| `pointing-check/` | 32 pointing tables; row order corresponds to M01-M19 |

The 74 reconstructed inputs total 42,719,344 bytes (42.7 MB). The VLA cube and mask use lossless gzip compression. Supporting inputs are bundled in `data-parts/supporting_inputs.zip`; `python prepare_data.py` unpacks and verifies them.

The large VLA H I cube is distributed separately as the **HCG100_HI.fits.gz** asset of the **v1.0.0 GitHub Release** (41,212,957 bytes). Place it in `VLA/` before running the notebooks. `python prepare_data.py --download` downloads and verifies it automatically.

**Release asset status:** Available in [v1.0.0](https://github.com/qz-yu/HCG100-sidelobe-correction/releases/tag/v1.0.0). The published asset size and SHA-256 digest match `data_downloads.json`.

Original uncompressed sizes and hashes and the packaged-file hashes are in `input_manifest.json`. The VLA input files decompress to exactly the original bytes. Beam and mask construction are upstream inputs to this release.

## Outputs and checks

Notebook 02 writes the four source maps and their sum to `sidelobe-result-fits/VLAgrid_method_FASTgrid_M0/`. Notebook 03 writes source and total sidelobe cubes and the corrected, unsmoothed cube to `sidelobe-result-fits/VLAgrid_method_channel_subtraction/`.

The final product is:

`HCG100-2fov-HIcube-pix1.4arcmin-20kms-1-mJyBeam-1-VLA-sidelobe-sub-smooth.fits`

After running all three notebooks, use `python validate.py --outputs` to check input hashes, pointing assignments, channel normalization, integrated model recovery, subtraction, and smoothing. See `VALIDATION.md` for the checks performed for this release and `CHANGES.md` for edits to the supplied notebooks.
