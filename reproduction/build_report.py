#!/usr/bin/env python3
"""Summarise actual notebook outputs without changing the official results."""
import hashlib,json,re,shutil
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results/reproduction'
RECEIPTS=ROOT/'reproduction/receipts'


def markdown_table(frame):
 def fmt(v):
  if isinstance(v,(float,np.floating)): return f'{v:.4f}'
  return str(v)
 return '\n'.join(['| '+' | '.join(map(str,frame.columns))+' |','| '+' | '.join(['---']*len(frame.columns))+' |']+['| '+' | '.join(fmt(v) for v in row)+' |' for row in frame.itertuples(index=False,name=None)])


def main():
 RECEIPTS.mkdir(parents=True,exist_ok=True)
 verification=json.loads((OUT/'verification.json').read_text())
 assert verification['status']=='passed'
 rows=[]; common=[]
 for panel,scope,cells in [('3a','Within slide',12719),('3c','Across slide',76779)]:
  d=pd.read_csv(OUT/f'Figure3/{panel}.csv')
  columns=[c for c in d if c.endswith(' PCC')]
  for n,label in [(20,'Top 20 SVGs'),(50,'Top 50 SVGs'),(len(d),'All shared genes')]:
   part=d.loc[d['SVG rank']<=n]
   same=part.loc[np.isfinite(part[columns]).all(axis=1)]
   for c in columns:
    rows.append({'Scope':scope,'Gene set':label,'Method':c[:-4],'Median PCC':float(part[c].median()),'Defined genes':int(part[c].notna().sum()),'Requested genes':len(part),'Cells':cells})
    common.append({'Scope':scope,'Gene set':label,'Method':c[:-4],'Median PCC':float(same[c].median()),'Common defined genes':len(same),'Requested genes':len(part),'Cells':cells})
 official=pd.DataFrame(rows); equal=pd.DataFrame(common)
 official.to_csv(OUT/'expression_summary.csv',index=False); equal.to_csv(OUT/'expression_common_defined_summary.csv',index=False)
 external=[]
 for panel in ['3h','3i']:
  d=pd.read_csv(OUT/f'Figure3/{panel}.csv')
  for n,label in [(20,'Top 20 SVGs'),(50,'Top 50 SVGs'),(len(d),'All shared genes')]:
   v=d.loc[d.SVG_rank<=n]
   external.append({'Dataset':d.Dataset.iloc[0],'Gene set':label,'Median PCC':float(v.PCC.median()),'Defined genes':int(v.PCC.notna().sum()),'Requested genes':len(v),'Matched cells':int(d.PCC_n.iloc[0])})
 external=pd.DataFrame(external); external.to_csv(OUT/'external_summary.csv',index=False)
 imputation=[]
 for metric,file in [('PCC','fig4a_pcc.csv'),('SSIM','fig4a_ssim.csv'),('RMSE','fig4a_rmse.csv')]:
  d=pd.read_csv(OUT/'Figure4'/file)
  for method in ['GHIST+','PASTA','SpatialEx','SpatialEx+']:
   v=d[f'{method} {metric}']; imputation.append({'Metric':metric,'Method':method,'Median':v.median(),'Defined slide-gene pairs':int(v.notna().sum())})
 imputation=pd.DataFrame(imputation); imputation.to_csv(OUT/'imputation_summary.csv',index=False)
 colors={'GHIST+':'#D1495B','GHIST':'#EDAE49','SpatialEx':'#4C8DB8','Phoenix':'#2A9D8F'}
 fig,axes=plt.subplots(1,3,figsize=(12,4.6),sharey=True)
 for ax,gene_set in zip(axes,['Top 20 SVGs','Top 50 SVGs','All shared genes']):
  table=official.loc[(official.Scope=='Across slide') & (official['Gene set']==gene_set)]
  ax.bar(table.Method,table['Median PCC'],color=[colors[m] for m in table.Method],width=.7)
  for i,(_,row) in enumerate(table.iterrows()):
   ax.text(i,row['Median PCC']+.018,f'{row["Median PCC"]:.3f}\nn={row["Defined genes"]}/{row["Requested genes"]}',ha='center',va='bottom',fontsize=8)
  ax.set_title(gene_set,fontsize=11); ax.tick_params(axis='x',rotation=30); ax.set_ylim(0,.85)
  ax.spines[['top','right']].set_visible(False)
 axes[0].set_ylabel('Median per-gene Pearson correlation')
 fig.suptitle('Across-slide evaluation of released predictions: 76,779 matched cells',fontsize=13)
 fig.text(.5,.015,'Official per-method summaries; undefined correlations excluded. See equal-gene comparison in the report.',ha='center',fontsize=8)
 fig.tight_layout(rect=(0,.05,1,.93)); fig.savefig(OUT/'across_slide_summary.png',dpi=180); fig.savefig(OUT/'across_slide_summary.pdf'); plt.close(fig)
 manifest=json.loads((ROOT/'reproduction/data_manifest.json').read_text())
 runs=json.loads((OUT/'run_manifest.json').read_text())
 coverage=pd.read_csv(OUT/'Figure5/coverage_summary.csv')
 small=official.loc[(official['Gene set']=='Top 50 SVGs')].drop(columns=['Cells'])
 common_small=equal.loc[(equal.Scope=='Across slide') & (equal['Gene set']=='Top 50 SVGs')].drop(columns=['Scope','Cells'])
 baseline=imputation.pivot(index='Method',columns='Metric',values='Median').reset_index()
 audit=json.loads((OUT/'Figure4/figure4f_audit.json').read_text())
 ready=json.loads((OUT/'model_readiness.json').read_text())
 report=f'''# GHIST+ reproduction report

Run completed **15 October 2026 (Asia/Kolkata)** in this workspace.

**Status: all four released figure notebooks executed successfully. This is a reproduction of evaluation from released predictions and summary tables. Fresh H&E-to-RNA inference and training have not been performed.**

## What was run

| Notebook | Scope reproduced | Evidence level |
| --- | --- | --- |
| Figure2 | Panels b–f: ROI/pathway maps and correlations | Replotting authors' released CSV summaries; original input checksums passed |
| Figure3 | Panels a–i, including ATERA and Breast5K | Recomputed metrics/maps from released cell-level predictions and measured RNA |
| Figure4 | Panels a, b, c, e and f provided by notebook | Recomputed imputation/structure metrics and VQ/composition ablation from released predictions |
| Figure5 | Three transfer directions | Replotting authors' paired per-gene PCC table; coverage checked; raw transfer predictions not rerun |

The original code cells were executed unchanged. One final cell per notebook preserves temporary plots and tables. No synthetic cells, RNA profiles or images were generated. The published randomized VQ control is an existing released ablation on real cells, not newly generated synthetic data.

## Source and environment

- Official source: https://github.com/SydneyBioX/GHIST_plus
- Source commit: `{runs[0]['upstream_commit']}`.
- Official data: https://huggingface.co/datasets/SydneyBioX/GHIST-Plus-bundle
- Dataset revision: `{manifest['revision']}`.
- {len(manifest['downloads'])} downloaded files; {sum(i['bytes'] for i in manifest['downloads'])/1e9:.3f} GB. Every file was checked against the source LFS SHA-256 or Git blob hash.
- Apple Silicon macOS, 16 GiB RAM; analysis environment uses Python 3.12.14. This differs from the authors' Python 3.10/CUDA training environment.
- Main numeric package versions follow the authors' tested versions. Matplotlib inline backend pinned to 0.1.7 for compatibility with Matplotlib 3.8.2.
- Full package versions: [pip_freeze.txt](pip_freeze.txt). Training seeds/configs remain in the unchanged upstream files.

## Expression prediction results

PCC measures correlation, **not percentage accuracy**. SVGs are genes ranked as spatially variable in measured target RNA. These selected genes do not represent the full gene panel.

{markdown_table(small)}

The within-slide cohort contains 12,719 matched cells and 280 shared genes. The across-slide cohort contains 76,779 cells and 276 shared genes. The numbers above reproduce the notebook's per-method median, excluding undefined correlations. In the across-slide top-50 comparison, GHIST has only 42 defined correlations; GHIST+, SpatialEx and Phoenix each have 50.

For a separate comparison on **exactly the same genes with finite PCC for every method**, the across-slide top-50 results are:

{markdown_table(common_small)}

This additional coverage diagnostic does not replace or modify the official notebook results. The all-shared-gene across-slide median for GHIST+ is **0.2031** (276/276 defined), well below its selected top-50 median of 0.5508. In the within-slide top-20 SVG group, original GHIST slightly exceeds GHIST+ (0.7746 vs 0.7623). Phoenix also has lower/better correlation-matrix distance in the released across-slide comparison. These results do not support a claim that GHIST+ wins every task.

![Across-slide released-prediction evaluation](across_slide_summary.png)

## External datasets

{markdown_table(external)}

ATERA has 168,335 unique matched cells and 413 shared genes; 403 genes have defined PCC. Breast5K has 564,723 matched cells and 258 shared genes, all with defined PCC. These are the released breast multi-slide model's external evaluations, not the Figure5 PanCancer comparisons. Counts refer to cells, not independent patients.

## Gene imputation and ablation

Median metrics over **45 shared held-out slide–gene cases across five real slides**:

{markdown_table(baseline)}

Higher PCC/SSIM is better; lower RMSE is better. These use the authors' saved prediction transforms and evaluation protocol. The 45 cases include different genes on different slides; they are not 45 independent patients.

The VQ/composition ablation used {audit['common_cells']:,} cells and {audit['common_genes']} genes. The full-model top-20 median PCC was 0.7586, versus 0.5730 when the VQ/composition pathway was off and 0.5810 for the published randomized control. Top-50 medians were 0.6425, 0.5183 and 0.5122. This evaluates saved frozen-checkpoint interventions; no control was retrained here.

The Figure3 within-slide cohort has 12,719 cells, while this Figure4 ablation uses 13,621. Their scores must not be treated as identical evaluation cohorts.

## PanCancer transfer: coverage matters

Figure5 is generated from a released per-gene summary table. Its paired comparisons have only **5/13** genes for pancreas top-20/top-50, **5/9** for breast, and **2/4** for lung. The union target panels are larger. Plotted pair counts should not be described as 20 or 50 jointly evaluable genes in every direction. Raw per-cell PanCancer predictions for these case studies were not independently re-evaluated by this notebook run.

{markdown_table(coverage[['Direction','Gene_set','target_n','paired_n','single_median_paired','pancancer_median_paired']])}

## Verification

- Every original source code cell matches the executed notebook exactly.
- All four notebooks completed with no execution errors; 16 displayed plot outputs were saved.
- The notebooks' cohort, gene-coverage and input assertions passed.
- Independent SciPy PCC checks recomputed three SVG-ranked genes across four methods from the real count/prediction files on 76,779 matched cells. Maximum difference from saved notebook PCC: {verification['maximum_absolute_pcc_difference']:.3g}.
- External cohort checks passed for ATERA and Breast5K; the five-slide imputation check passed.
- [verification.json](verification.json) contains the check receipts; [run_manifest.json](run_manifest.json) contains timings and source hashes. This verifies execution/numeric consistency, not independent model training or exact correspondence to every value in a manuscript version.

## Remaining requirements for a full model reproduction

The official figure bundle does not contain the aligned training H&E images, nuclei masks and matched-nuclei files required by the training configuration. Preflight found {len(ready['missing_input_fields'])} missing path fields for `config_multi_breast.json` (some fields are related inputs for the same slide).

UNI2-H is a gated encoder. The GHIST+ released checkpoints exclude its tensors; its official weights require an approved Hugging Face account and institutional-email access: https://huggingface.co/MahmoodLab/UNI2-h . No encoder weights were substituted and no access application was submitted. The current analysis environment has no PyTorch/CUDA installation; this Mac has no NVIDIA CUDA GPU. The authors' tested training setup uses Linux/CUDA.

Consequently, **new image predictions, released-checkpoint inference, fresh training, baseline retraining, and patient-level replication remain unfinished**. The original training/inference code and configurations are cloned locally, and the readiness check names each missing input. A full rerun must retain the published slide splits and settings; a reduced Mac run would be a different experiment.

## Open the results

- [Figure2 executed notebook](Figure2/Figure2.executed.ipynb)
- [Figure3 executed notebook](Figure3/Figure3.executed.ipynb)
- [Figure4 executed notebook](Figure4/Figure4.executed.ipynb)
- [Figure5 executed notebook](Figure5/Figure5.executed.ipynb)
- [Figure3 across-slide plot](Figure3/3c.png)
- [Figure4 ablation plot](Figure4/figure4f_vq_ablation.png)
- [Figure5 transfer plot](Figure5/Figure5.png)

The original Figure4 ablation strip titles overlap slightly in the authors' plotting layout; its numerical data are intact. The separate summary figure above has a readable layout without altering the original notebook.
'''
 (OUT/'RELEASED_NOTEBOOK_REVIEW.md').write_text(report)
 (ROOT/'reproduction/archive/released_notebook_review_regenerated.md').write_text(re.sub(r'\]\((?!https?://)([^)]+)\)', r'](../../results/reproduction/\1)',report))
 shutil.copy2(OUT/'across_slide_summary.png',RECEIPTS/'across_slide_summary.png')
 for name in ['expression_summary.csv','expression_common_defined_summary.csv','external_summary.csv','imputation_summary.csv','verification.json','run_manifest.json','environment.json','model_readiness.json','pip_freeze.txt']:
  shutil.copy2(OUT/name,RECEIPTS/name)
 for fig in [2,3,4,5]: shutil.copy2(OUT/f'Figure{fig}/execution.json',RECEIPTS/f'Figure{fig}_execution.json')
 freeze=(OUT/'pip_freeze.txt').read_text(); (ROOT/'reproduction/requirements-analysis.lock.txt').write_text('# Exact macOS ARM64 / Python 3.12 analysis environment; CUDA training uses a separate environment.\n'+freeze)
 artifact_hashes={str(p.relative_to(OUT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.rglob('*') if p.is_file() and not any(part.startswith('.') for part in p.relative_to(OUT).parts) and p.name!='artifact_checksums.json'}
 (OUT/'artifact_checksums.json').write_text(json.dumps(artifact_hashes,indent=2)+'\n')
 shutil.copy2(OUT/'artifact_checksums.json',RECEIPTS/'artifact_checksums.json')
 print('Report and compact receipts saved:',OUT/'RELEASED_NOTEBOOK_REVIEW.md')

if __name__=='__main__': main()
