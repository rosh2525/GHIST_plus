#!/usr/bin/env python3
"""Build a manuscript-scoped report; never substitute paper targets for measured results."""
import hashlib,json,shutil
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
SPEC=ROOT/'reproduction/paper_protocol'
OLD=ROOT/'results/reproduction'
OUT=ROOT/'results/paper_only'


def table(df):
 def cell(v):
  if pd.isna(v): return 'Not reproduced'
  if isinstance(v,float): return f'{v:.4f}'
  return str(v).replace('|','/')
 return '\n'.join(['| '+' | '.join(df.columns)+' |','| '+' | '.join(['---']*len(df.columns))+' |']+['| '+' | '.join(cell(v) for v in row)+' |' for row in df.itertuples(index=False,name=None)])


def main():
 protocol=json.loads((SPEC/'protocol.json').read_text()); paper=protocol['paper']
 pdf=ROOT/paper['local_pdf']
 assert pdf.is_file() and hashlib.sha256(pdf.read_bytes()).hexdigest()==paper['pdf_sha256'],'Run reproduction/fetch_paper.py first'
 values=pd.read_csv(SPEC/'paper_values.csv')
 allowed={m['name'] for m in protocol['metrics']}
 assert set(values.metric)<=allowed
 OUT.mkdir(parents=True,exist_ok=True)
 values.to_csv(OUT/'paper_reported_values.csv',index=False)
 # Preserve provenance of the completed released-code evaluation.
 for f in ['Figure2','Figure3','Figure4']:
  receipt=json.loads((OLD/f/'execution.json').read_text())
  assert receipt['status']=='passed'
  assert hashlib.sha256((ROOT/(f+'.ipynb')).read_bytes()).hexdigest()==receipt['source_notebook_sha256']
  shutil.copy2(OLD/f/'execution.json',OUT/(f+'_released_execution.json'))
 within=pd.read_csv(OLD/'Figure3/3a.csv'); across=pd.read_csv(OLD/'Figure3/3c.csv')
 rows=[]
 for label,frame,n in [('Within slide',within,20),('Within slide',within,50),('Across slide',across,20),('Across slide',across,50)]:
  for method in ['GHIST+','GHIST','SpatialEx']:
   v=frame.loc[frame['SVG rank']<=n,method+' PCC']
   rows.append({'Experiment':label,'Gene set':f'Top {n} SVGs','Method':method,'Metric':'PCC','Released-data value':float(v.median()),'Defined genes':int(v.notna().sum()),'Requested genes':n,'Status':'Released-code result; manuscript equivalence unverified'})
 released=pd.DataFrame(rows); released.to_csv(OUT/'released_code_paper_metric_results.csv',index=False)
 # Manuscript panel names, not the later repository's renumbered panels.
 for name,source in [('Figure3a_PCC_per_gene.csv','Figure3/3a.csv'),('Figure3b_SSIM_per_gene.csv','Figure3/3b.csv'),('Figure3c_PCC_per_gene.csv','Figure3/3c.csv'),('Figure3d_CMD.csv','Figure3/3f.csv')]:
  d=pd.read_csv(OLD/source)
  d=d.drop(columns=[c for c in d if c.startswith('Phoenix ')])
  if 'Gene subset' in d: d=d.loc[d['Gene subset'].isin(['All genes','Top 20 SVGs','Top 50 SVGs','Top 20 non-SVGs','Top 50 non-SVGs'])]
  d.to_csv(OUT/name,index=False)
 auc=pd.read_csv(OLD/'Figure3/3d.csv');auc.loc[auc.Model.isin(['GHIST+','GHIST','SpatialEx'])].to_csv(OUT/'Figure3e_AUROC.csv',index=False)
 iou=pd.read_csv(OLD/'Figure3/3e.csv');iou.loc[iou.Source.isin(['GHIST+','GHIST','SpatialEx']) & iou.Class.ne('Unassigned')].to_csv(OUT/'Figure3f_Spatial_IoU.csv',index=False)
 # Imputation values stay outside the strict manuscript results: scaling/cohort differ.
 normalization={
  'paper_equation':'Eq.15: log(1 + counts/library_size * 10000)',
  'released_notebook':'Figure4 load_imputation_slide uses log1p(raw counts), GHIST+ predictions / 2 and log1p(PASTA).',
  'consequence':'The released RMSE cannot be presented as a reproduction of manuscript RMSE. Missing full measured panels/library sizes must be recovered first.',
  'heldout_cohort':'Manuscript masks top ten genes per slide. Current released evaluation has 45 common slide-gene cases over five slides; original eligibility and fold manifests are unresolved.'}
 (OUT/'normalization_audit.json').write_text(json.dumps(normalization,indent=2)+'\n')
 comparisons=[]
 for key,method,expected in [('within_ghistplus','GHIST+',.744),('within_ghist','GHIST',.714),('within_spatialex','SpatialEx',.521)]:
  actual=float(within.loc[within['SVG rank']<=20,method+' PCC'].median())
  comparisons.append({'Experiment':'Within slide, top 20 SVGs','Method':method,'Paper-reported PCC':expected,'Released-code PCC':actual,'Status':'Does not reproduce manuscript value'})
 comparison=pd.DataFrame(comparisons);comparison.to_csv(OUT/'paper_vs_released.csv',index=False)
 panel_map=[
  ['2b-f','Figure2.ipynb','Paper-relevant plots exist; numerical mismatch prevents an exact reproduction claim.'],
  ['3a,b,c','Figure3 tables a,b,c','Current released cohort; three paper methods retained; five-fold paper result not reconstructed.'],
  ['3d,e,f','Figure3 tables f,d,e respectively','CMD/AUROC/IoU remapped to manuscript numbering; results remain released-code results.'],
  ['4a,b,d','Figure4.ipynb','Available partial outputs; paper normalization, eligibility and complete baselines unresolved.'],
  ['4c','Later Figure4 table f','Current frozen-checkpoint intervention has not been established as the manuscript no-VQ experiment.'],
  ['5a','Current Figure5.ipynb is different','Paper-specific PancreasADC and PancreasADC2 transfer has not been reproduced.'],
  ['5b-f','Figure5_new.Rmd','GTEx sample objects, compartments and molecular predictions are missing; analyses not run.']
 ]
 pd.DataFrame(panel_map,columns=['Paper panels','Current repository counterpart','Reproduction status']).to_csv(OUT/'figure_mapping.csv',index=False)
 metrics=pd.DataFrame([{'Metric':m['name'],'Paper page':','.join(map(str,m['paper_pages'])),'Meaning':m['unit']} for m in protocol['metrics']])
 cell_labels=values.loc[values.id.str.startswith('figure3'),['cohort','metric','value']].pivot(index='cohort',columns='metric',values='value').reset_index().rename(columns={'cohort':'Cell type'})
 selected=values.loc[values.id.isin(['within_ghistplus','within_ghist','within_spatialex','across_ghistplus','imputation_pcc_1','imputation_pcc_2','imputation_rmse_1','imputation_rmse_2','imputation_tissue','pancreas2_single20','pancreas2_multi20','pancreas2_single50','pancreas2_multi50','pancreas1_single20','pancreas1_multi20']),['experiment','method','metric','cohort','value','pdf_page']]
 report=f'''# GHIST+ manuscript-only reproduction

Scope updated **22 October 2026**. Source: *{paper['title']}*, bioRxiv v1, posted 12 June 2026, DOI {paper['doi']}.

**Exact numerical reproduction of the manuscript has not been achieved.** The completed local run reproduces the authors' later released notebooks. Those notebooks differ from the manuscript in reported values, panel numbering, normalization and configuration. A successful notebook run does not establish a successful paper reproduction.

## Metrics defined in this paper

{table(metrics)}

SSIM and AUROC should not be called percent accuracy. TISSUE is an error/variability score, not confidence-interval coverage. The paper also presents downstream cell proportions, ECM scores, regression slopes and P-values for GTEx; these are not single-cell RNA prediction accuracy.

## Values actually reported by the paper

The following are manuscript targets, **not new measurements from our machine**. Values come from the Results text on PDF pages 5–6. The full transcription, including pathway results, is in `paper_reported_values.csv`.

{table(selected)}

The manuscript lists two PCC and two RMSE imputation values together. Their complete protocol assignments are not specified in that sentence; they are retained as separate reported results rather than assigning an invented average or cohort. SSIM, CMD, AUROC, recall and IoU are included in the metric specification; no precise scalar was invented for a boxplot without a reported numerical value.

## Cell-type values printed on the paper figure

These GHIST+ values are rounded labels from Figure3e/f (PDF page30), not results obtained by training here. The corresponding file retains their figure provenance.

{table(cell_labels)}

## Current released-code results differ

{table(comparison)}

The current released within-slide cohort has 12,719 matched cells. The paper describes five 80/20 spatial folds, while the checked-in config supplies only one 0–20% validation division. These differences have not been reconciled, so the released values cannot be labelled as matched manuscript results. We do not tune or replace observed results just to obtain a target number.

## Figure correspondence

{table(pd.DataFrame(panel_map,columns=['Paper panels','Current counterpart','Status']))}

The main paper compares GHIST+, GHIST and SpatialEx for expression prediction. Phoenix, ATERA/Breast5K external benchmarks, the extra common-finite-gene analysis and the SpatialEx+ arm are excluded from this manuscript-only report. No new architecture or improvement experiment is introduced. Earlier broader output is archived separately.

## Specific issues preventing a full match

1. **Normalization.** Manuscript RMSE uses library-size-normalized log expression (Eq.15). The released imputation notebook uses log1p of raw counts. Its previous RMSE is therefore excluded from the manuscript result table.
2. **Folds and settings.** Paper: five folds, learning rate 0.001, weight decay 0.0001. Released config: one division, learning rate 0.00115, weight decay 0.00001. Templates restore the two stated optimizer settings; they do not invent the missing original fold assignments.
3. **Panel selection.** The imputation release evaluates 45 common slide–gene cases; the paper describes masking the top ten genes on each of five slides. Original exclusions and complete measured panels are needed before claiming the same cohort.
4. **Variants.** The later VQ/composition intervention and manuscript ablation are not established as the same experiment. Baseline outputs for SpaGE, SpaIM, stImpute and Tangram described in the text are also absent locally.
5. **Figure5.** The manuscript includes pancreas transfer and 77-sample normal-breast/age analyses. The current Python notebook's three transfer directions do not reproduce that full figure. The R script refers to unavailable author-side objects.
6. **Model execution.** Aligned H&E/masks/matched inputs, approved UNI2-H weights and a suitable model environment are still needed. No fresh model training or image inference has run.

The paper itself has caption inconsistencies: Figure3b plots HVGs but its caption says non-SVGs; Figure4c plots non-SVGs but its caption says HVGs. These are recorded rather than silently equating HVGs with non-SVGs. The Methods define a signed TISSUE error while the released notebook aggregates absolute errors; that convention also needs reconciliation.

## Workflow

Read [WORKFLOW.md](WORKFLOW.md) for the exact data → preprocessing → model → experiments → metrics sequence and commands. The original PDF is preserved unmodified at `{paper['local_pdf']}` with SHA-256 `{paper['pdf_sha256']}`.

Execution receipts and source hashes are in `paper_scope_receipt.json`. All exported numeric rows use metrics mentioned in this manuscript; every source is labelled as manuscript text or later released-code evaluation.
'''
 (ROOT/'reproduction/REPORT.md').write_text(report)
 (OUT/'REPORT.md').write_text(report.replace('(WORKFLOW.md)','(../../reproduction/WORKFLOW.md)'))
 shutil.copy2(OUT/'paper_vs_released.csv',SPEC/'paper_vs_released.csv')
 receipt={'paper':paper,'scope':'manuscript-only','paper_value_rows':len(values),'allowed_metrics':sorted(allowed),'fresh_training_run':False,'fresh_image_inference_run':False,'exact_paper_reproduction_achieved':False,'notebook_evaluation_completed_previously':True,'current_sources':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [OLD/'Figure3/3a.csv',OLD/'Figure3/3c.csv',OLD/'Figure3/3f.csv',OLD/'Figure3/3d.csv',OLD/'Figure3/3e.csv']},'outputs_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in OUT.iterdir() if p.is_file() and p.name!='paper_scope_receipt.json'}}
 (OUT/'paper_scope_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 (SPEC/'paper_scope_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('Paper-only report generated; exact manuscript reproduction remains incomplete:',OUT/'REPORT.md')

if __name__=='__main__':main()
