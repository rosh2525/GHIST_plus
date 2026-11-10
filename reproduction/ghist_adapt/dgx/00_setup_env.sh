#!/usr/bin/env bash
# One-time environment on the DGX (Linux + CUDA). Mirrors the authors' tested stack (README: torch 2.6.0 / CUDA 12.4,
# Python 3.10) plus what the ghist_adapt scripts need.
set -euo pipefail
ENV=${ENV:-ghist_adapt}
if ! command -v conda >/dev/null; then echo "conda/mamba not found: load your cluster's module first (e.g. module load anaconda)"; exit 1; fi
conda create -y -n "$ENV" python=3.10
eval "$(conda shell.bash hook)"; conda activate "$ENV"
pip install torch==2.6.0 torchvision==0.21.0 --index-url https://download.pytorch.org/whl/cu124
pip install numpy==1.26.4 pandas==2.3.2 scipy==1.15.3 scikit-learn==1.7.1 matplotlib==3.8.2 tabulate tqdm natsort h5py \
            tifffile imagecodecs "zarr<3" imageio opencv-python-headless pillow timm==0.9.12 "huggingface_hub>=0.34" torchstain
python - <<'PY'
import torch; print("torch", torch.__version__, "cuda", torch.cuda.is_available(), torch.cuda.device_count(), "GPU(s)")
PY
echo "Next: huggingface-cli login (token with access to MahmoodLab/UNI2-h), then dgx/01_fetch.sh"
