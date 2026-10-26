# Changes from the supplied notebooks

The three supplied notebooks are the sole source for this release. Their source-file SHA-256 hashes are recorded in `input_manifest.json`.

- Removed saved outputs, execution counts, empty cells and local notebook metadata; added short workflow introductions.
- Renamed `01_prepare_vla_scaled_truth_maps.ipynb.ipynb` to `01_prepare_vla_scaled_truth_maps.ipynb`.
- Replaced `WCS.proj_plane_pixel_scales() * u.deg` with `u.Quantity(WCS.proj_plane_pixel_scales(), u.deg)` where needed. The WCS method returns angular quantities in current Astropy; multiplying by degrees again produces squared angular units.
- Corrected the header comment describing the target convolution beam: it comes from the FAST header.
- Corrected notebook 02's misleading recentering docstring. Its existing implementation uses the stored kernel arrays without a CRPIX shift.
- Corrected notebook 02 output `BUNIT` to `mJy/beam km/s` instead of inheriting the input surface-brightness unit. The output maps retain the reference geometry and beam keywords from the VLA header; those beam keywords do not describe an effective Gaussian resolution for the sidelobe model.
- Corrected notebook 03's output comment to identify the smoothed VLA aperture spectra.
- Made the final smoothing loop use `nchan` instead of a hard-coded 140-channel velocity array. The supplied cube has 140 channels, so this does not change the processed channels.

Source positions, aperture radii, fixed scale factors, velocity selections, beam resampling and normalization, convolution, pointing averaging, spectral weighting, subtraction and smoothing parameters are retained.

- Packaged the VLA cube and mask as lossless .fits.gz files and updated the two notebook input paths. Decompressed bytes match the original inputs exactly.
- Bundled the supporting inputs in a verified ZIP archive. The large compressed VLA cube is provided separately as a GitHub Release asset; `prepare_data.py` unpacks the archive and verifies the cube.
