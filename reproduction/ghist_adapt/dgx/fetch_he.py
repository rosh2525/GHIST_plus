"""Download a slide's 10x H&E OME-TIFF and alignment CSV listed in dgx/he_sources.json.
Usage: python dgx/fetch_he.py breast5k   ->  $GHIST_BUNDLE_ROOT/he/<slide>/{he.ome.tif, he_alignment.csv, SHA256SUMS}"""
import hashlib
import json
import os
import sys
import urllib.request
from pathlib import Path

src = json.loads((Path(__file__).parent / "he_sources.json").read_text())
slide = sys.argv[1]
s = src[slide]
if not s["he_ome_tif_url"] or not s["alignment_csv_url"]:
    sys.exit(f"Fill he_ome_tif_url and alignment_csv_url for '{slide}' in dgx/he_sources.json (see {s['page']})")
out = Path(os.environ["GHIST_BUNDLE_ROOT"]) / "he" / slide
out.mkdir(parents=True, exist_ok=True)
sums = []
for url, name in ((s["he_ome_tif_url"], "he.ome.tif"), (s["alignment_csv_url"], "he_alignment.csv")):
    dst = out / name
    if not dst.exists():
        print("downloading", url)
        urllib.request.urlretrieve(url, dst)
    h = hashlib.sha256(dst.read_bytes()).hexdigest()
    sums.append(f"{h}  {name}  {url}")
(out / "SHA256SUMS").write_text("\n".join(sums) + "\n")
print("\n".join(sums))
