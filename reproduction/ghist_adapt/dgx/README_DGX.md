# DGX runbook (P4)

Everything below runs from `reproduction/ghist_adapt/`. The protocols are frozen JSON files in the parent folder.
Order matters: the cheap, certain jobs come first.

## 0. One-time setup (about 30 min)
```bash
git clone <your copy of this repo> && cd GHIST_plus/reproduction/ghist_adapt
bash dgx/00_setup_env.sh                       # conda env "ghist_adapt": torch 2.6 / CUDA 12.4 (authors' stack)
conda activate ghist_adapt
huggingface-cli login                          # token whose account accepted MahmoodLab/UNI2-h terms (only needed for P4b)
export GHIST_BUNDLE_ROOT=/scratch/$USER/ghist_bundle GHIST_ADAPT_CACHE=/scratch/$USER/ghist_adapt_cache
export GHIST_V3_FETCH=$GHIST_BUNDLE_ROOT
bash dgx/01_fetch.sh                           # about 20 GB; WITH_CHECKPOINTS=1 also fetches the released checkpoints
```

## 1. P4a: second-checkpoint confirmation (CPU, about 2–4 h, no GPU or UNI2-h needed). Run this first.
The pan-cancer GHIST+ checkpoint predicted the **same Breast5K slide** (its held-out slide 11) with its own,
larger gene vocabulary. It tests whether the zero-shot soft-expectation fix, HAC and P2 hold for a different model.
Protocols are unchanged (config_v4/config_v5).
```bash
sbatch dgx/02_confirm_pancancer.sh             # or: bash dgx/02_confirm_pancancer.sh
```
Outputs go to `outputs_v4/breast5k_pan/` and `outputs_v5/breast5k_pan/`, plus the logs `p4a_p1.log` and `p4a_p2.log`.

## 2. P4b: image-aware correction (GPU)
Needs (a) UNI2-h access and (b) the 10x **post-Xenium H&E OME-TIFF + H&E alignment CSV** for the slide. The 10x
website blocked automated fetching from here (HTTP 429). Copy the two links from the dataset page's download
section into `dgx/he_sources.json`:
- Breast5K: https://www.10xgenomics.com/datasets/xenium-prime-ffpe-human-breast-cancer
- BreastCancer2 (Rep 2): https://www.10xgenomics.com/products/xenium-in-situ/preview-dataset-human-breast

```bash
python dgx/fetch_he.py breast5k                # downloads to $GHIST_BUNDLE_ROOT/he/<slide>/ from he_sources.json
# QA first: the green dots must sit on nuclei; otherwise flip --align-direction or --cell-units
python dgx/extract_uni2h.py --he $GHIST_BUNDLE_ROOT/he/breast5k/he.ome.tif --align $GHIST_BUNDLE_ROOT/he/breast5k/he_alignment.csv \
       --cells $GHIST_BUNDLE_ROOT/evaluation_data/breast5k/cell_coords_histology_filtered.csv --cell-units um \
       --out /scratch/$USER/emb/breast5k_uni2h --qa-only
# full run (about 565k cells; roughly 30–60 min on one A100/H100)
python dgx/extract_uni2h.py ... --out /scratch/$USER/emb/breast5k_uni2h --batch 512 --workers 16
# freeze the protocol: set "frozen" in config_v6.json (date/time) BEFORE this:
python dgx/p4_run.py breast5k --emb /scratch/$USER/emb/breast5k_uni2h
```
Caveat: `cell_coords_histology_filtered.csv` is in micrometres. Whether it is in the Xenium frame or already in the
H&E frame is exactly what the QA overlay checks. Do not run the full extraction until the overlay is right.

## 3. Not in this package (and why)
- **Re-running GHIST+ inference, or the authors' `train_tma_select.py` head-to-head.** Both need each slide's
  `aligned_he_image.tif`, HoverNet nucleus mask and matched-nuclei table. The authors did not release these. Rebuilding
  them means replicating the GHIST preprocessing (raw Xenium output, H&E registration, HoverNet on GPU, nucleus–cell
  matching), roughly a week. It is only worth it if P4a/P4b results justify it. BreastCancer2 hurdle outputs need
  the same pipeline.

## Resource summary
| Job | Node | Time (estimate) | Disk |
|---|---|---|---|
| fetch | any | 20–40 min | about 20 GB |
| P4a | CPU, 16 cores, 128 GB RAM | 2–4 h | about 15 GB cache |
| UNI2-h extraction | 1 GPU (A100/H100), 16 CPU workers | 0.5–1 h per slide | H&E 2–10 GB, embeddings about 3.5 GB per slide (float16) |
| P4b run | CPU, 64 GB | about 30 min per slide | small |
