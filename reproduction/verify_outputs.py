#!/usr/bin/env python3
"""Verify notebook provenance/cohorts and independently recompute selected PCCs."""
import hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import pearsonr
ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'results'/'reproduction'
BUNDLE=ROOT/'local_artifacts'/'bundle'

def keyed_csv(file,genes):
 columns=list(pd.read_csv(file,nrows=0).columns)
 key=next((k for k in ['c_id','cell_id','cellID','Cell ID','id','V1'] if k in columns),columns[0])
 frame=pd.read_csv(file,usecols=[key,*genes]).set_index(key)
 frame.index=frame.index.astype(np.int64)
 assert frame.index.is_unique, str(file)
 return frame

def main():
 checks=[]
 for fig in [2,3,4,5]:
  folder=OUT/f'Figure{fig}'; metadata=json.loads((folder/'execution.json').read_text())
  original=ROOT/f'Figure{fig}.ipynb'; executed=json.loads((folder/f'Figure{fig}.executed.ipynb').read_text())
  source=json.loads(original.read_text())
  assert metadata['status']=='passed'
  assert hashlib.sha256(original.read_bytes()).hexdigest()==metadata['source_notebook_sha256']
  assert len(executed['cells'])==len(source['cells'])+1
  for before,after in zip(source['cells'],executed['cells']):
   assert ''.join(before['source'])==''.join(after['source']),f'Original code changed: Figure{fig}'
   if after['cell_type']=='code': assert after['execution_count'] is not None
   assert not any(out['output_type']=='error' for out in after.get('outputs',[]))
  checks.append({'check':f'Figure{fig} original code, execution counts and no errors','passed':True})
 # Verify real-data alignment and PCC independently of notebook helper functions.
 pcc_table=pd.read_csv(OUT/'Figure3'/'3c.csv')
 chosen=pcc_table.sort_values('SVG rank').iloc[[0,19,49]]['Gene'].tolist()
 base=BUNDLE/'other_models'/'across_slide'
 paths={'GHIST+':BUNDLE/'GHIST_plus/predictions/breast_multi/test_slide3_pred_expr_scaled.csv',
        'GHIST':base/'GHIST_pred_expr.csv','SpatialEx':base/'SpatialEx_UNI2H_pred_expr.csv','Phoenix':base/'Phoenix_pred_expr.csv'}
 actual=keyed_csv(BUNDLE/'evaluation_data/breast2/cell_gene_matrix_filtered.csv',chosen)
 coord=keyed_csv(BUNDLE/'evaluation_data/breast2/cell_coords.csv',[])
 labels=keyed_csv(BUNDLE/'evaluation_data/breast2/cell_type_filtered.csv',[])
 pred={name:keyed_csv(file,chosen) for name,file in paths.items()}
 ids=actual.index.intersection(coord.index).intersection(labels.index)
 for frame in pred.values(): ids=ids.intersection(frame.index)
 assert len(ids)==76779
 differences=[]
 for name,frame in pred.items():
  for gene in chosen:
   measured=pearsonr(np.log1p(actual.loc[ids,gene]),frame.loc[ids,gene]).statistic
   saved=float(pcc_table.set_index('Gene').loc[gene,name+' PCC'])
   assert np.isclose(measured,saved,atol=1e-12,rtol=1e-12,equal_nan=True)
   differences.append({'method':name,'gene':gene,'recalculated_pcc':float(measured) if np.isfinite(measured) else None,'saved_pcc':saved if np.isfinite(saved) else None,'absolute_difference':float(abs(measured-saved)) if np.isfinite(measured) else None})
 checks.append({'check':'Independent SciPy PCC on 76,779 matched cells, three SVG-ranked genes and four methods','passed':True})
 for name,rows,valid,cells in [('3h.csv',413,403,168335),('3i.csv',258,258,564723)]:
  table=pd.read_csv(OUT/'Figure3'/name)
  assert len(table)==rows and table.PCC.notna().sum()==valid
  assert table.PCC_n.eq(cells).all()
  checks.append({'check':f'{name}: expected genes, defined PCCs and unique-cell cohort','passed':True})
 imputation=pd.read_csv(OUT/'Figure4/fig4a_pcc.csv')
 assert len(imputation)==45 and imputation.Slide.nunique()==5
 assert not imputation.duplicated(['Slide','Gene']).any()
 checks.append({'check':'45 distinct held-out slide-gene pairs over five slides','passed':True})
 receipt={'status':'passed','checks':checks,'independent_pcc_checks':differences,'maximum_absolute_pcc_difference':max(i['absolute_difference'] for i in differences if i['absolute_difference'] is not None)}
 (OUT/'verification.json').write_text(json.dumps(receipt,indent=2,allow_nan=False)+'\n')
 print(json.dumps(receipt,indent=2,allow_nan=False))

if __name__=='__main__': main()
