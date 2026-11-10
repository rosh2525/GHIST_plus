"""Per-cell pathology-foundation-model embeddings from a 10x post-Xenium H&E OME-TIFF (Apple Metal, CUDA or CPU).

  python extract_features.py qa  --slide breast5k        # choose + verify the coordinate transform (no model)
  python extract_features.py run --slide breast5k --encoder phikon

Encoders: 'phikon' (owkin/phikon, ViT-B/16, open download, Owkin non-commercial licence) or 'uni2h' (MahmoodLab/UNI2-h
via the authors' model/backbone.py, gated). Crop: 224 px at the pyramid level nearest 0.5 um/px (resampled to exactly
0.5 um/px = 112 um field of view, the magnification both encoders were trained at), centred on each cell.
Output: <cache>/<slide>_<encoder>.npy (float16, n_cells x 2*dim = [CLS, mean patch token]) + _ids.npy.

Coordinates: cell (x, y) in um -> Xenium pixels (/ 0.2125) -> H&E level-0 pixels through the 10x alignment matrix.
The qa step scores every candidate transform (matrix or its inverse, x/y as given or swapped) by how much darker
(more haematoxylin) the H&E is at the cell centroids than 30 um away, picks the best, and writes an overlay PNG.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
REPO = HERE.parents[1]
HE = REPO / "local_artifacts" / "he"
from slides import CACHE, SLIDES  # noqa: E402

XENIUM_PX_UM = 0.2125
TARGET_UM = 0.5


def open_levels(path):
    import tifffile
    tf = tifffile.TiffFile(path)
    s = tf.series[0]
    px = None
    try:
        import re
        m = re.search(r'PhysicalSizeX="([0-9.eE+-]+)"', tf.ome_metadata or "")
        px = float(m.group(1)) if m else None
    except Exception:
        pass
    levels = []
    for i, lv in enumerate(s.levels):
        shp = lv.shape
        chw = len(shp) == 3 and shp[0] in (3, 4) and shp[-1] not in (3, 4)
        H, W = (shp[1], shp[2]) if chw else (shp[0], shp[1])
        levels.append(dict(i=i, H=H, W=W, chw=chw))
    for lv in levels:
        lv["down"] = levels[0]["W"] / lv["W"]
    return tf, levels, px


def level_array(path, i):
    import tifffile, zarr
    z = zarr.open(tifffile.imread(path, aszarr=True), mode="r")
    return z[str(i)] if hasattr(z, "keys") and str(i) in z else z


def cell_table(slide):
    c = pd.read_csv(SLIDES[slide]["coords"])
    c = c.rename(columns={c.columns[0]: "id", "x_coord": "x", "y_coord": "y"})
    return c["id"].to_numpy(np.int64), c[["x", "y"]].to_numpy(np.float64) * SLIDES[slide]["um_per_unit"]


def transform(xy_um, A, variant):
    inv, swap = variant
    xy = xy_um[:, ::-1] if swap else xy_um
    xp = xy / XENIUM_PX_UM
    M = np.linalg.inv(A) if inv else A
    h = np.c_[xp, np.ones(len(xp))] @ M.T
    return h[:, :2] / h[:, 2:3]


def read_patch(arr, chw, x0, y0, w, h, H, W):
    out = np.full((h, w, 3), 255, np.uint8)
    xs, ys, xe, ye = max(x0, 0), max(y0, 0), min(x0 + w, W), min(y0 + h, H)
    if xe > xs and ye > ys:
        t = np.asarray(arr[:3, ys:ye, xs:xe]) if chw else np.asarray(arr[ys:ye, xs:xe, :3])
        out[ys - y0:ye - y0, xs - x0:xe - x0] = np.moveaxis(t, 0, -1) if chw else t
    return out


def qa(slide):
    d = HE / slide
    tf, levels, px0 = open_levels(d / "he.ome.tif")
    A = np.loadtxt(d / "he_alignment.csv", delimiter=",")
    ids, xy = cell_table(slide)
    rng = np.random.default_rng(0)
    pick = rng.choice(len(xy), 3000, replace=False)
    # work on a ~2 um/px level for speed
    lv = min(levels, key=lambda l: abs(l["down"] * (px0 or 0.27) - 2.0))
    arr = level_array(d / "he.ome.tif", lv["i"])
    img_small = None
    scores = {}
    off = 30.0 / ((px0 or 0.27) * lv["down"])            # 30 um in level pixels
    for variant in [(True, False), (False, False), (True, True), (False, True)]:
        p = transform(xy[pick], A, variant) / lv["down"]
        ok = (p[:, 0] > off) & (p[:, 1] > off) & (p[:, 0] < lv["W"] - off) & (p[:, 1] < lv["H"] - off)
        if ok.mean() < 0.5:
            scores[str(variant)] = {"inside": float(ok.mean()), "contrast": None}
            continue
        vals, ctrl = [], []
        for (x, y) in p[ok][:600]:
            patch = read_patch(arr, lv["chw"], int(x) - 1, int(y) - 1, 3, 3, lv["H"], lv["W"]).mean()
            far = read_patch(arr, lv["chw"], int(x + off) - 1, int(y) - 1, 3, 3, lv["H"], lv["W"]).mean()
            vals.append(patch); ctrl.append(far)
        scores[str(variant)] = {"inside": float(ok.mean()), "contrast": float(np.mean(ctrl) - np.mean(vals)),
                                "brightness_at_cells": float(np.mean(vals))}
        print(variant, scores[str(variant)], flush=True)
    best = max((k for k in scores if scores[k]["contrast"] is not None), key=lambda k: scores[k]["contrast"])
    variant = tuple(v == "True" for v in best.strip("()").replace(" ", "").split(","))
    info = {"slide": slide, "level0_um_per_px": px0, "levels": [(l["i"], l["W"], l["H"]) for l in levels],
            "scores": scores, "chosen": {"invert_matrix": variant[0], "swap_xy": variant[1]}}
    # overlay at full resolution for 6 cells
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from scipy.spatial import cKDTree
    arr0 = level_array(d / "he.ome.tif", 0)
    P0 = transform(xy, A, variant)
    tree = cKDTree(P0)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for ax, i in zip(axes.flat, rng.choice(len(xy), 6, replace=False)):
        w = 300
        img = read_patch(arr0, levels[0]["chw"], int(P0[i, 0]) - w // 2, int(P0[i, 1]) - w // 2, w, w, levels[0]["H"], levels[0]["W"])
        nb = tree.query_ball_point(P0[i], w / 2)
        ax.imshow(img); ax.scatter(P0[nb, 0] - P0[i, 0] + w / 2, P0[nb, 1] - P0[i, 1] + w / 2, s=10, c="lime")
        ax.set_title(f"cell {ids[i]}"); ax.axis("off")
    fig.suptitle(f"{slide}: centroids (green) on H&E, transform {info['chosen']}"); fig.tight_layout()
    fig.savefig(d / "qa_overlay.png", dpi=100)
    (d / "qa.json").write_text(json.dumps(info, indent=2))
    print(json.dumps(info["chosen"]), "->", d / "qa_overlay.png")


def load_encoder(name, dev):
    import torch
    if name == "phikon":
        from transformers import ViTModel
        m = ViTModel.from_pretrained("owkin/phikon", add_pooling_layer=False)
        fwd = lambda x: m(pixel_values=x).last_hidden_state                      # [CLS] + 196 patches
        n_skip = 1
    elif name == "uni2h":
        import sys; sys.path.insert(0, str(REPO))
        from model.backbone import Uni2HEncoder
        m = Uni2HEncoder(pretrained=True)
        fwd = lambda x: (lambda t: t["x"] if isinstance(t, dict) else t)(m.vit.forward_features(x))
        n_skip = 1 + 8
    else:
        raise ValueError(name)
    dtype = torch.float16 if dev.type in ("mps", "cuda") else torch.float32
    m = m.to(dev, dtype=dtype).eval()
    return fwd, n_skip, dtype


def run(slide, encoder, batch, strip_rows):
    import torch
    d = HE / slide
    info = json.loads((d / "qa.json").read_text())
    variant = (info["chosen"]["invert_matrix"], info["chosen"]["swap_xy"])
    tf, levels, px0 = open_levels(d / "he.ome.tif")
    px0 = px0 or info["level0_um_per_px"] or 0.2738
    lv = min((l for l in levels if l["down"] * px0 <= TARGET_UM * 1.05), key=lambda l: abs(l["down"] * px0 - TARGET_UM))
    lv_um = lv["down"] * px0
    win = int(round(224 * TARGET_UM / lv_um))                  # crop at this level that covers 112 um
    arr = level_array(d / "he.ome.tif", lv["i"])
    A = np.loadtxt(d / "he_alignment.csv", delimiter=",")
    ids, xy = cell_table(slide)
    P = transform(xy, A, variant) / lv["down"]
    order = np.argsort(P[:, 1], kind="stable")
    dev = torch.device("cuda" if torch.cuda.is_available() else "mps" if torch.backends.mps.is_available() else "cpu")
    fwd, n_skip, dtype = load_encoder(encoder, dev)
    mean = torch.tensor([0.485, 0.456, 0.406], device=dev, dtype=dtype).view(1, 3, 1, 1)
    std = torch.tensor([0.229, 0.224, 0.225], device=dev, dtype=dtype).view(1, 3, 1, 1)
    out = None
    t0 = time.time(); done = 0
    print(f"{slide}: level {lv['i']} ({lv_um:.3f} um/px), crop {win} px -> 224, {len(ids)} cells, device {dev}", flush=True)
    i = 0
    while i < len(order):
        y0 = int(P[order[i], 1]) - win // 2
        j = i
        while j < len(order) and int(P[order[j], 1]) - win // 2 < y0 + strip_rows:
            j += 1
        sel = order[i:j]
        ys = max(y0, 0); ye = min(int(P[sel, 1].max()) + win // 2 + 1, lv["H"])
        strip = np.asarray(arr[:3, ys:ye, :]) if lv["chw"] else np.asarray(arr[ys:ye, :, :3])
        if lv["chw"]:
            strip = np.moveaxis(strip, 0, -1)
        for b in range(0, len(sel), batch):
            cb = sel[b:b + batch]
            crops = np.full((len(cb), win, win, 3), 255, np.uint8)
            for k, c in enumerate(cb):
                x0 = int(round(P[c, 0])) - win // 2; yy = int(round(P[c, 1])) - win // 2 - ys
                xs, xe = max(x0, 0), min(x0 + win, lv["W"]); y_s, y_e = max(yy, 0), min(yy + win, strip.shape[0])
                if xe > xs and y_e > y_s:
                    crops[k, y_s - yy:y_e - yy, xs - x0:xe - x0] = strip[y_s:y_e, xs:xe]
            x = torch.from_numpy(crops).to(dev).permute(0, 3, 1, 2).to(dtype).div_(255)
            if win != 224:
                x = torch.nn.functional.interpolate(x, size=224, mode="bilinear", antialias=True, align_corners=False)
            with torch.inference_mode():
                tok = fwd((x - mean) / std)
                emb = torch.cat([tok[:, 0], tok[:, n_skip:].mean(1)], 1).float().cpu().numpy().astype(np.float16)
            if out is None:
                out = np.zeros((len(ids), emb.shape[1]), np.float16)
            out[cb] = emb
            done += len(cb)
        i = j
        el = time.time() - t0
        print(f"  {done}/{len(ids)} cells, {done/el:.0f}/s, eta {(len(ids)-done)/max(done/el,1e-9)/3600:.1f} h", flush=True)
    dst = CACHE / f"{slide}_{encoder}"
    np.save(str(dst) + ".npy", out); np.save(str(dst) + "_ids.npy", ids)
    json.dump({"slide": slide, "encoder": encoder, "level": lv["i"], "um_per_px": lv_um, "crop_px": win,
               "fov_um": 224 * TARGET_UM, "transform": info["chosen"], "seconds": round(time.time() - t0)},
              open(str(dst) + ".json", "w"), indent=2)
    print("saved", dst, out.shape)


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("step", choices=["qa", "run"])
    ap.add_argument("--slide", default="breast5k")
    ap.add_argument("--encoder", default="phikon", choices=["phikon", "uni2h"])
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--strip-rows", type=int, default=2048)
    a = ap.parse_args()
    qa(a.slide) if a.step == "qa" else run(a.slide, a.encoder, a.batch, a.strip_rows)
