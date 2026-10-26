# Release validation

Validated on 2026-10-26 using Python 3.12 on Windows. Numerical package versions are recorded in `requirements-validated.txt`.

All code cells in notebooks 01, 02 and 03 were executed in order, with a fresh Python namespace for each notebook and non-interactive plotting. The published notebooks have their saved outputs cleared.

`python validate.py --outputs` passed:

- SHA-256 and byte counts for all 74 input files.
- Syntax and notebook-format checks for all three notebooks (notebook-format validation was run separately).
- 32 pointing tables with 19 beam positions each: 608 hits on 560 FAST pixels.
- FAST cube shape: 140 × 34 × 28, with 20 km/s channels.
- All four source spectra: `sum(channel_weight × dv) = 1`.
- Per-source channel models agree with the integrated maps times their independently reconstructed spectral weights.
- Integrated channel models recover the original source sidelobe maps. The largest absolute recovery difference is 6.15 × 10⁻⁷ mJy beam⁻¹ km s⁻¹, consistent with float32 storage.
- The saved total cube equals the sum of source cubes.
- The corrected cube equals the input cube minus the total sidelobe model.
- Every smoothed channel agrees with the specified convolution; the input NaN mask is retained and the output beam is labelled 4′ × 4′.

Machine-readable results are in `validation_results.json`. These checks verify execution and internal numerical consistency of the supplied method. They do not independently establish the accuracy of the input beam models, source templates or adopted scale factors.

The VLA mask carries the non-standard FITS unit label `source_ID`, which spectral-cube warns about. It is used only as the Boolean selection `mask > 0`; the warning did not prevent execution.
