"""P4b: per-cell UNI2-h embeddings from a 10x post-Xenium H&E image (GPU).

Uses the authors' encoder (model/backbone.py::Uni2HEncoder): the same weights and architecture GHIST+ uses.
Needs Hugging Face access to MahmoodLab/UNI2-h (huggingface-cli login).

Steps:
  1) --qa-only first. It writes a PNG of 6 random cells with their centroids drawn on the H&E crop. CHECK IT: the
     dots must sit on nuclei. If they don't, fix --align-direction / --cell-units before extracting.
  2) the full run writes <out>.npy (float16, n_cells x 3072 = [CLS, mean patch token]) and <out>_ids.npy.

Coordinates: cells (id, x, y) are converted to Xenium-image pixels (x_um / xenium_px_um), then to H&E pixels using
the 10x H&E alignment matrix (3x3 affine). 10x's alignment CSV maps H&E pixels to Xenium pixels by default
('he_to_xenium'); the script inverts it. Verify with --qa-only.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np
import pandas as pd

REPO = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO))


def he_coords(cells: pd.DataFrame, align: np.ndarray, cell_units: str, xenium_px_um: float, direction: str):
    xy = cells[["x", "y"]].to_numpy(np.float64)
    if cell_units == "um":
        xy = xy / xenium_px_um                      # -> Xenium pixels
    if cell_units != "he_px":
        M = np.linalg.inv(align) if direction == "he_to_xenium" else align
        h = np.c_[xy, np.ones(len(xy))] @ M.T
        xy = h[:, :2] / h[:, 2:3]
    return xy                                        # H&E level-0 pixels (x, y)


def open_he(path, level):
    import tifffile, zarr
    store = tifffile.imread(path, aszarr=True)
    z = zarr.open(store, mode="r")
    arr = z[str(level)] if isinstance(z, zarr.hierarchy.Group) else z
    return arr                                       # (Y, X, C) or (C, Y, X)


def crop(arr, cx, cy, size):
    chw = arr.shape[0] in (3, 4) and arr.shape[-1] not in (3, 4)
    H, W = (arr.shape[1], arr.shape[2]) if chw else (arr.shape[0], arr.shape[1])
    x0, y0 = int(round(cx)) - size // 2, int(round(cy)) - size // 2
    xs, ys = max(x0, 0), max(y0, 0); xe, ye = min(x0 + size, W), min(y0 + size, H)
    out = np.full((size, size, 3), 255, dtype=np.uint8)
    if xe > xs and ye > ys:
        tile = arr[:3, ys:ye, xs:xe] if chw else arr[ys:ye, xs:xe, :3]
        tile = np.asarray(tile); tile = np.moveaxis(tile, 0, -1) if chw else tile
        out[ys - y0:ye - y0, xs - x0:xe - x0] = tile[..., :3]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--he", required=True, help="post-Xenium H&E OME-TIFF")
    ap.add_argument("--cells", required=True, help="CSV with cell_id,x_coord,y_coord (e.g. evaluation_data/.../cell_coords*.csv)")
    ap.add_argument("--align", help="10x H&E alignment CSV (3x3 matrix); omit when --cell-units he_px")
    ap.add_argument("--align-direction", choices=["he_to_xenium", "xenium_to_he"], default="he_to_xenium")
    ap.add_argument("--cell-units", choices=["um", "xenium_px", "he_px"], default="um")
    ap.add_argument("--xenium-px-um", type=float, default=0.2125)
    ap.add_argument("--level", type=int, default=0)
    ap.add_argument("--patch", type=int, default=224, help="crop size in H&E pixels (fed to UNI2-h at 224)")
    ap.add_argument("--batch", type=int, default=256)
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--out", required=True)
    ap.add_argument("--qa-only", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    args = ap.parse_args()

    c = pd.read_csv(args.cells)
    c = c.rename(columns={c.columns[0]: "id", "x_coord": "x", "y_coord": "y"})
    align = np.loadtxt(args.align, delimiter=",") if args.align else np.eye(3)
    xy = he_coords(c, align, args.cell_units, args.xenium_px_um, args.align_direction)
    arr = open_he(args.he, args.level)
    print(f"H&E level {args.level} shape {arr.shape}; cell H&E px x[{xy[:,0].min():.0f},{xy[:,0].max():.0f}] y[{xy[:,1].min():.0f},{xy[:,1].max():.0f}]")

    # QA overlay
    import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
    from scipy.spatial import cKDTree
    rng = np.random.default_rng(0); tree = cKDTree(xy)
    fig, axes = plt.subplots(2, 3, figsize=(12, 8))
    for ax, i in zip(axes.flat, rng.choice(len(xy), 6, replace=False)):
        img = crop(arr, xy[i, 0], xy[i, 1], args.patch)
        nb = tree.query_ball_point(xy[i], args.patch / 2)
        ax.imshow(img); ax.scatter(xy[nb, 0] - xy[i, 0] + args.patch / 2, xy[nb, 1] - xy[i, 1] + args.patch / 2, s=8, c="lime")
        ax.set_title(f"cell {c.id.iloc[i]}"); ax.axis("off")
    fig.suptitle("QA: green dots must sit on nuclei"); fig.tight_layout(); fig.savefig(args.out + "_qa.png", dpi=110)
    print("QA overlay:", args.out + "_qa.png")
    if args.qa_only:
        return

    import torch
    from torch.utils.data import DataLoader, Dataset
    from model.backbone import Uni2HEncoder

    n = len(xy) if not args.limit else min(args.limit, len(xy))
    mean = torch.tensor([0.485, 0.456, 0.406]).view(1, 3, 1, 1); std = torch.tensor([0.229, 0.224, 0.225]).view(1, 3, 1, 1)

    class Crops(Dataset):
        arr = None                                    # opened lazily, once per DataLoader worker

        def __len__(self): return n

        def __getitem__(self, i):
            if self.arr is None:
                self.arr = open_he(args.he, args.level)
            return torch.from_numpy(crop(self.arr, xy[i, 0], xy[i, 1], args.patch))

    dev = torch.device("cuda")
    enc = Uni2HEncoder(pretrained=True).to(dev).eval()
    out = np.zeros((n, 2 * enc.out_channels), dtype=np.float16)
    dl = DataLoader(Crops(), batch_size=args.batch, num_workers=args.workers, pin_memory=True)
    k = 0
    with torch.inference_mode(), torch.autocast("cuda", dtype=torch.bfloat16):
        for b in dl:
            x = b.permute(0, 3, 1, 2).float().div(255)
            if args.patch != 224:
                x = torch.nn.functional.interpolate(x, size=224, mode="bilinear", antialias=True)
            x = ((x - mean) / std).to(dev, non_blocking=True)
            tok = enc.vit.forward_features(x)
            tok = tok["x"] if isinstance(tok, dict) else tok
            cls, patches = tok[:, 0], tok[:, 1 + 8:]           # [CLS] + 8 register tokens + patch tokens
            out[k:k + len(x)] = torch.cat([cls, patches.mean(1)], 1).float().cpu().numpy().astype(np.float16)
            k += len(x)
            if (k // args.batch) % 50 == 0:
                print(f"{k}/{n}", flush=True)
    np.save(args.out + ".npy", out); np.save(args.out + "_ids.npy", c.id.to_numpy()[:n])
    print("saved", args.out + ".npy", out.shape)


if __name__ == "__main__":
    main()
