"""Unpack supporting inputs and verify or download the VLA H I cube."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import zipfile

ROOT = Path(__file__).resolve().parent

def prepare_data(download=False):
    manifest = json.loads((ROOT / "data_downloads.json").read_text())
    expected = {f["path"]: f for f in json.loads((ROOT / "input_manifest.json").read_text())["files"]}
    support = manifest["support_archive"]
    archive = ROOT / support["path"]
    if hashlib.sha256(archive.read_bytes()).hexdigest() != support["sha256"]:
        raise ValueError("Supporting input archive failed its checksum.")
    with zipfile.ZipFile(archive) as bundle:
        for name in bundle.namelist():
            destination = ROOT / name
            if name not in expected or not destination.resolve().is_relative_to(ROOT.resolve()):
                raise ValueError("Unexpected archive member: " + name)
            data = bundle.read(name)
            if hashlib.sha256(data).hexdigest() != expected[name]["sha256"]:
                raise ValueError("Invalid input: " + name)
            destination.parent.mkdir(parents=True, exist_ok=True)
            destination.write_bytes(data)
    print("Supporting input files unpacked and verified.")
    cube = manifest["vla_cube"]
    destination = ROOT / cube["path"]
    if not destination.exists():
        if not download:
            raise SystemExit("Place HCG100_HI.fits.gz in VLA/, or run: python prepare_data.py --download\nRelease asset: " + cube["url"])
        temporary = destination.with_suffix(".gz.downloading")
        try:
            with urllib.request.urlopen(cube["url"], timeout=60) as response, temporary.open("wb") as handle:
                while block := response.read(1024*1024):
                    handle.write(block)
            if temporary.stat().st_size != cube["size_bytes"] or hashlib.sha256(temporary.read_bytes()).hexdigest() != cube["sha256"]:
                raise ValueError("Downloaded cube failed its checksum.")
            temporary.replace(destination)
        except Exception:
            temporary.unlink(missing_ok=True)
            raise
    if destination.stat().st_size != cube["size_bytes"] or hashlib.sha256(destination.read_bytes()).hexdigest() != cube["sha256"]:
        raise ValueError("VLA H I cube failed its checksum.")
    print("VLA H I cube verified.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--download", action="store_true", help="download the VLA cube from the v1.0.0 GitHub Release if missing")
    prepare_data(download=parser.parse_args().download)
