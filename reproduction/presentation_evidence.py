#!/usr/bin/env python3
"""Print concise evidence from real files for presentation terminal captures."""
import argparse,json,subprocess,sys
from pathlib import Path
import pandas as pd
ROOT=Path(__file__).resolve().parents[1]
OLD=ROOT/'results/reproduction'
BUNDLE=ROOT/'local_artifacts/bundle'


def load(name):return json.loads((ROOT/name).read_text())

def ids_genes(p):
 cols=list(pd.read_csv(p,nrows=0).columns)
 key=next((k for k in ['c_id','cell_id','cellID','Cell ID','id','V1'] if k in cols),cols[0])
 ids=set(pd.read_csv(p,usecols=[key])[key].astype('int64'))
 return ids,set(cols)-{key,'x','y','x_coord','y_coord','x coordinate','y coordinate'}

def dataset_summary():
 gt=BUNDLE/'evaluation_data/breast2/cell_gene_matrix_filtered.csv'
 coord=BUNDLE/'evaluation_data/breast2/cell_coords.csv'
 labels=BUNDLE/'evaluation_data/breast2/cell_type_filtered.csv'
 gtids,gtgenes=ids_genes(gt);coids,_=ids_genes(coord);laids,_=ids_genes(labels)
 within=[BUNDLE/'GHIST_plus/predictions/breast_single/trainval_slide3_pred_expr_scaled.csv',BUNDLE/'other_models/within_slide/GHIST_pred_expr.csv',BUNDLE/'other_models/within_slide/SpatialEx_UNI2H_pred_expr.csv']
 across=[BUNDLE/'GHIST_plus/predictions/breast_multi/test_slide3_pred_expr_scaled.csv',BUNDLE/'other_models/across_slide/GHIST_pred_expr.csv',BUNDLE/'other_models/across_slide/SpatialEx_UNI2H_pred_expr.csv',BUNDLE/'other_models/across_slide/Phoenix_pred_expr.csv']
 def cohort(files,include_labels):
  cells=gtids&coids; genes=gtgenes.copy()
  if include_labels:cells&=laids
  for f in files:
   ci,ge=ids_genes(f);cells&=ci;genes&=ge
  return len(cells),len(genes)
 wi,wg=cohort(within,False);ai,ag=cohort(across,True)
 imp=pd.read_csv(OLD/'Figure4/fig4b_cmd.csv')
 result={'BreastCancer2_input_cells':len(gtids),'BreastCancer2_genes':len(gtgenes),'within_cells':wi,'within_genes':wg,'across_cells':ai,'across_genes':ag,'imputation_slides':len(imp),'imputation_slide_gene_cases':int(imp['Holdout gene count'].sum()),'imputation_per_slide':imp[['Slide','Cell count','Holdout gene count']].to_dict('records')}
 print(f'BreastCancer2 input    {len(gtids):,} cells   {len(gtgenes)} genes')
 print(f'Within-slide cohort   {wi:,} cells   {wg} shared genes')
 print(f'Across-slide cohort   {ai:,} cells   {ag} shared genes')
 print(f'Imputation cohort     {len(imp)} slides       {int(imp["Holdout gene count"].sum())} held-out slide-gene cases')
 print('Source: official processed 10x Xenium data + released predictions')
 (ROOT/'results/paper_only/presentation_dataset_summary.json').write_text(json.dumps(result,indent=2)+'\n')


def checks():
 result=subprocess.run([sys.executable,str(ROOT/'reproduction/verify_outputs.py')],cwd=ROOT,capture_output=True,text=True)
 if result.returncode:print(result.stdout+result.stderr);raise SystemExit(result.returncode)
 manifest=load('reproduction/data_manifest.json');receipt=load('results/reproduction/verification.json')
 runs=load('results/reproduction/run_manifest.json')
 print(f'Recorded source-verified downloads: {len(manifest["downloads"])} files, {sum(i["bytes"] for i in manifest["downloads"])/1e9:.3f} GB')
 print('Released notebook execution, 15 Oct 2026:')
 ordered=sorted(runs,key=lambda r:r['figure'])
 for pair in [ordered[:2],ordered[2:]]:print('  '+'     '.join(f'Figure{r["figure"]}: {r["status"].upper()}' for r in pair))
 print(f'Fresh provenance / cell-ID / PCC checks: {receipt["status"].upper()}')
 print('Fresh model training: false   Fresh image inference: false')


def compare():
 d=pd.read_csv(ROOT/'results/paper_only/paper_vs_released.csv')
 for _,row in d.iterrows():print(f'{row["Method"]:12} Paper PCC: {row["Paper-reported PCC"]:.3f}   Local PCC: {row["Released-code PCC"]:.4f}')


def status():
 result=subprocess.run([sys.executable,str(ROOT/'reproduction/preflight.py')],cwd=ROOT,capture_output=True,text=True)
 if result.returncode not in [0,2]:print(result.stdout+result.stderr);raise SystemExit(result.returncode)
 r=load('results/reproduction/model_readiness.json')
 cfg=load('configs/config_multi_breast.json');spec=load('reproduction/paper_protocol/protocol.json')
 print(f'CUDA available: {r["cuda_available"]}    UNI2-H weights supplied: {r["uni2_weights_supplied"]}')
 print(f'Missing aligned-data path fields: {len(r["missing_input_fields"])}')
 print(f'Validation divisions in release: {len(cfg["regions_val"]["divisions"])}    Folds in paper: {spec["training"]["within_slide_folds"]}')
 print('Prepared: manuscript settings + experiment templates')
 print('Pending: original folds, model inputs, encoder, fresh runs')
 print('Exact manuscript match remains incomplete')

if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('mode',choices=['data','checks','compare','status']);a=p.parse_args();globals()[{'data':'dataset_summary'}.get(a.mode,a.mode)]()
