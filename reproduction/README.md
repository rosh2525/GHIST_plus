# Local GHIST+ reproduction

**All four official analysis notebooks ran successfully on 15 October 2026. Fresh model training/image inference has not been performed.**

Open [the reproduction report](../results/reproduction/REPRODUCTION_REPORT.md) for outcomes, denominators, limitations and file links. Compact numerical receipts are versioned in [receipts](receipts/); downloaded data and full generated notebooks/plots are kept locally under ignored directories.

## Pinned sources

- Code: SydneyBioX/GHIST_plus, commit `c5c895174c19903cc27d2661934b9b726bdb9306`.
- Data: SydneyBioX/GHIST-Plus-bundle, revision `37a937677c5a8ebe5c7b8fbd486fbc5bb50ce1de`.
- Full source and checksum records: [data_manifest.json](data_manifest.json).

The upstream model, configurations and notebook code were preserved. The runner only appends output persistence code. No synthetic data is used. Figure2 and Figure5 replot published summaries; Figure3 and Figure4 recompute metrics from released real-cell predictions/ground truth.

## Repeat the analysis

Run from the `GHIST_plus` repository root. This setup uses a Python 3.12 analysis environment, independently of the authors' Python 3.10/CUDA model environment.

```sh
# Create analysis environment. Set GHIST_REPRO_PYTHON to a Python 3.12 executable if needed.
sh reproduction/setup_analysis.sh

# Official files only; about 10.9 GB; downloads are resumed and source hashes checked.
../work/ghist_plus_reproduction/.venv/bin/python reproduction/fetch_bundle.py --scope all

# Run official Figures 2, 4, 5, 3 serially to limit memory.
../work/ghist_plus_reproduction/.venv/bin/python reproduction/run_figures.py --figures 2 4 5 3

# Independent provenance/cohort/PCC verification.
../work/ghist_plus_reproduction/.venv/bin/python reproduction/verify_outputs.py

# This exits 2 until real image inputs, encoder weights and CUDA are available.
../work/ghist_plus_reproduction/.venv/bin/python reproduction/preflight.py

# Build the report after the preflight receipt has been written.
MPLCONFIGDIR=results/reproduction/.matplotlib ../work/ghist_plus_reproduction/.venv/bin/python reproduction/build_report.py
```

For exact analysis versions, install `requirements-analysis.lock.txt` in a Python 3.12 macOS ARM64 environment. Other platforms should use `requirements-analysis.txt` and record their own lock. The analysis kernel specification/caches are stored under the results directory; no global Jupyter kernel is registered.

`fetch_bundle.py --scope core` plus `run_figures.py --core-only` explicitly omits external Figure3 panels h/i, for smaller machines. The completed run here used `--scope all` and **did not omit any official notebook cells**.

## Fresh model execution

The original [README](../README.md), [tutorial](../tutorial.ipynb), [train.py](../train.py), [inference tool](../tools/inference.py) and [configs](../configs/) are included. Follow the authors' CUDA/Linux setup in a separate environment.

Required inputs: approved UNI2-H encoder weights; a reconstructed released GHIST+ checkpoint or a new training run; aligned H&E images, nuclei masks, matched cell metadata, and measured RNA/labels for evaluation/training. The public figure bundle alone does not include the image preparation inputs. Run `preflight.py --data-root /actual/prepared/data --uni2-weights /actual/pytorch_model.bin` to list missing configuration paths.

The current Mac has 16 GiB RAM and no NVIDIA CUDA GPU. Neither checkpoint inference nor training was replaced by an untrained/random encoder or claimed from the released predictions. No cloud machines were provisioned and no paid jobs were started.

## Git and licensing

Local branch: `reproduction/ghist-plus-local`; `origin` remains the authors' repository. No changes were pushed upstream. Code retains the upstream GPL-3.0 license; processed data, comparison predictions and UNI2-H have their own upstream terms (see downloaded `local_artifacts/bundle/THIRD_PARTY_NOTICES.md`).
