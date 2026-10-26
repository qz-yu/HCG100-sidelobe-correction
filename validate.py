"""Validate the released inputs and, optionally, the generated correction."""
import argparse
import ast
import hashlib
import gzip
import json
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent


def check_inputs():
    manifest = json.loads((ROOT / "input_manifest.json").read_text())
    for item in manifest["files"]:
        path = ROOT / item["path"]
        assert path.stat().st_size == item["size_bytes"], item["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == item["sha256"], item["path"]
        if "uncompressed_sha256" in item:
            raw = gzip.decompress(path.read_bytes())
            assert len(raw) == item["uncompressed_size_bytes"]
            assert hashlib.sha256(raw).hexdigest() == item["uncompressed_sha256"]
    for path in sorted(ROOT.glob("0[123]_*.ipynb")):
        nb = json.loads(path.read_text(encoding="utf-8"))
        for i, cell in enumerate(nb["cells"]):
            if cell["cell_type"] == "code":
                ast.parse("".join(cell["source"]), filename=f"{path.name}:cell{i}")
    for i in range(1, 33):
        points = np.loadtxt(ROOT / f"pointing-check/HCG100_p{i}_radec.list")
        assert points.shape == (19, 2)
        assert np.isfinite(points).all()
    return {"input_files": len(manifest["files"]), "input_bytes": sum(x["size_bytes"] for x in manifest["files"]),
            "pointings": 32, "beams_per_pointing": 19, "input_hashes": "PASS", "notebook_syntax": "PASS"}


def check_outputs():
    import astropy.units as u
    from astropy.io import fits
    from astropy.wcs import WCS
    from astropy.coordinates import SkyCoord
    from astropy.convolution import convolve, Gaussian2DKernel
    from spectral_cube import SpectralCube
    from photutils.aperture import CircularAperture, aperture_photometry

    u.add_enabled_units(u.def_unit(["mJy/Beam"], represents=u.Jy/u.beam/1000))
    raw = SpectralCube.read(ROOT / "rawcube/HCG100-2fov-HIcube-pix1.4arcmin-20kms-1-mJyBeam-1.fits")
    data = raw.filled_data[:].value
    vel = raw.spectral_axis.to_value(u.km/u.s)
    dv = float(np.median(np.abs(np.diff(vel))))
    assert np.isclose(dv, 20.0)
    vla = SpectralCube.read(ROOT / "VLA/vla_smooth_3arcmin_on_fastgrid_5000_5800_mJyBeam.fits")
    vdata = np.nan_to_num(vla.filled_data[:].value, nan=0.0)
    vv = vla.spectral_axis.to_value(u.km/u.s)
    assert vdata.shape[1:] == data.shape[1:]
    sources = None
    nb = json.loads((ROOT / "01_prepare_vla_scaled_truth_maps.ipynb").read_text(encoding="utf-8"))
    for cell in nb["cells"]:
        if cell["cell_type"] == "code":
            for node in ast.parse("".join(cell["source"])).body:
                if isinstance(node, ast.Assign) and any(isinstance(t, ast.Name) and t.id == "sources" for t in node.targets):
                    sources = ast.literal_eval(node.value)
    assert sources is not None

    mdir = ROOT / "sidelobe-result-fits/VLAgrid_method_FASTgrid_M0"
    cdir = ROOT / "sidelobe-result-fits/VLAgrid_method_channel_subtraction"
    wcs = WCS(raw.header).celestial
    expected_weight = np.zeros(data.shape[1:])
    for i in range(1, 33):
        points = np.loadtxt(ROOT / f"pointing-check/HCG100_p{i}_radec.list")
        x, y = wcs.world_to_pixel(SkyCoord(points[:,0]*u.deg, points[:,1]*u.deg))
        for xi, yi in zip(np.round(x).astype(int), np.round(y).astype(int)):
            if 0 <= yi < data.shape[1] and 0 <= xi < data.shape[2]:
                expected_weight[yi, xi] += 1
    total_model = np.zeros_like(data, dtype=float)
    source_results = {}
    for name, src in sources.items():
        smap = fits.getdata(mdir / f"sidelobe_{name}_M0_FASTgrid.fits")
        model = fits.getdata(cdir / f"sidelobe_{name}_cube_scaled_to_channels.fits")
        weight = fits.getdata(mdir / f"weight_{name}_FASTgrid.fits")
        np.testing.assert_array_equal(weight, expected_weight)
        assert (smap[weight == 0] == 0).all()
        xy = wcs.world_to_pixel(SkyCoord(src["ra"]*u.deg, src["dec"]*u.deg))
        scales = u.Quantity(wcs.proj_plane_pixel_scales(), u.deg).to_value(u.arcmin)
        radius = src["radius_arcmin"] / np.sqrt(abs(scales[0]*scales[1]))
        aperture = CircularAperture(xy, r=radius)
        spec = np.array([float(aperture_photometry(img, aperture)["aperture_sum"][0])
                         if src["vmin"] <= v <= src["vmax"] else 0.0
                         for img, v in zip(vdata, vv)])
        spec = np.maximum(spec, 0)
        norm = spec.sum()*dv
        assert norm > 0
        weights = spec/norm
        full = np.zeros_like(vel)
        for i, vf in enumerate(vel):
            j = np.argmin(abs(vv-vf))
            if abs(vv[j]-vf) < 0.5*dv:
                full[i] = weights[j]
        np.testing.assert_allclose(full.sum()*dv, 1.0, rtol=1e-10)
        np.testing.assert_allclose(model, smap[None,:,:]*full[:,None,None], rtol=5e-6, atol=1e-5)
        recovered = model.sum(axis=0, dtype=float)*dv
        np.testing.assert_allclose(recovered, smap, rtol=5e-6, atol=1e-4)
        total_model += model
        source_results[name] = {"channel_weight_integral": float(full.sum()*dv),
                               "max_M0_recovery_error_mJy_beam_kms": float(np.max(abs(recovered-smap)))}
    saved_total = fits.getdata(cdir / "sidelobe_total_cube_scaled_to_channels.fits")
    np.testing.assert_allclose(saved_total, total_model, rtol=5e-6, atol=1e-5)
    corrected = fits.getdata(cdir / "HCG100_FAST_cube_VLAgrid_sidelobe_sub.fits")
    np.testing.assert_allclose(corrected, data-saved_total, rtol=5e-6, atol=1e-5, equal_nan=True)
    final_path = ROOT / "HCG100-2fov-HIcube-pix1.4arcmin-20kms-1-mJyBeam-1-VLA-sidelobe-sub-smooth.fits"
    final = fits.getdata(final_path)
    h = fits.getheader(final_path)
    assert np.isclose(h["BMAJ"]*60, 4.0) and np.isclose(h["BMIN"]*60, 4.0)
    assert np.array_equal(np.isnan(final), np.isnan(data))
    kernel = Gaussian2DKernel(x_stddev=0.8344, y_stddev=0.9631)
    for i in range(len(data)):
        expected = convolve(corrected[i], kernel, boundary="fill", fill_value=np.nan,
                            mask=np.isnan(data[i]), preserve_nan=True)
        np.testing.assert_allclose(final[i], expected, rtol=5e-6, atol=1e-5, equal_nan=True)
    return {"cube_shape": list(data.shape), "dv_kms": dv, "sampled_pixels": int((expected_weight>0).sum()),
            "pointing_hits_on_grid": int(expected_weight.sum()), "source_checks": source_results,
            "model_sum": "PASS", "subtraction": "PASS", "smoothing_and_mask": "PASS"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outputs", action="store_true", help="also validate products after running all notebooks")
    parser.add_argument("--json-output", type=Path, help="save the check results as JSON")
    args = parser.parse_args()
    result = {"inputs": check_inputs()}
    if args.outputs:
        result["outputs"] = check_outputs()
    text = json.dumps(result, indent=2)
    print(text)
    if args.json_output:
        args.json_output.write_text(text+"\n", encoding="utf-8")
