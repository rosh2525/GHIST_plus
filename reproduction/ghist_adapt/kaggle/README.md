# Running the remaining work on free Kaggle compute

Kaggle notebooks are free with a verified account: about 30 GB RAM and 4 CPU cores, with T4/P100 GPUs for about 30 h per week. That is enough for everything the DGX package was meant to do.

1. **Build the code bundle on the Mac:** `bash reproduction/ghist_adapt/kaggle/build_zip.sh` → `local_artifacts/ghist_adapt_kaggle.zip` (about 110 KB).
2. **On kaggle.com:** go to Datasets → New Dataset, upload the zip, name it `ghist-adapt-code`, and keep it private.
3. **Notebook:** Code → New Notebook → File → Import Notebook → `kaggle/ghist_adapt_kaggle.ipynb`. Add the dataset with *Add Data*. Settings: Internet on; Accelerator *GPU T4 x1* (Part B only).
4. **Part A (CPU, about 2–3 h):** run P1–P3 on the **pan-cancer checkpoint's** Breast5K predictions, a second released GHIST+ model, under the unchanged protocols. Download `ghist_adapt_results_P4a.zip` from the Output tab and unzip it into `reproduction/ghist_adapt/` on the Mac.
5. **Part B (GPU):**
   - Store your Hugging Face token as a Kaggle Secret named `HF_TOKEN`. The account must have accepted the MahmoodLab/UNI2-h terms; this usually needs an institutional email.
   - Paste the 10x H&E and alignment links into `dgx/he_sources.json` before zipping.
   - Check the QA image before running the full extraction and P4b.

The notebook copies results to `/kaggle/working`, so they survive the session. Kaggle stops idle interactive sessions; for the long UNI2-h extraction use *Save Version → Save & Run All* so it runs in the background (up to 12 h per run).
