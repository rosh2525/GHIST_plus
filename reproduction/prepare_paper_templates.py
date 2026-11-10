#!/usr/bin/env python3
"""Apply only explicit manuscript training settings, without inventing fold boundaries."""
import copy,json
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
DEST=ROOT/'reproduction/paper_protocol/config_templates'

def main():
 DEST.mkdir(parents=True,exist_ok=True)
 sources={
 'within_slide':'configs/breast2_ablation9/full.json',
 'breast_transfer':'configs/config_multi_breast.json',
 'imputation':'configs/config_gene_mask_imputer.json',
 'pancancer':'configs/config_pancancer.json'}
 receipt=[]
 for name,source in sources.items():
  d=json.loads((ROOT/source).read_text()); before=copy.deepcopy(d)
  d['training'].update(learning_rate=.001,weight_decay=.0001,total_epochs=30,batch_size=64,beta1=.9,beta2=.999)
  if name=='pancancer':d['regions_test']={'divisions':[[0.,1.]]}
  d['experiment_dirs']={'load_dir':'results/paper_model_runs/'+name,'model_dir':'models','metrics_dir':'metrics'}
  (DEST/(name+'.json')).write_text(json.dumps(d,indent=2)+'\n')
  receipt.append({'template':name+'.json','source':source,'optimizer_before':before['training'],'manuscript_override':{'learning_rate':.001,'weight_decay':.0001,'total_epochs':30,'batch_size':64,'beta1':.9,'beta2':.999},'validation_divisions_preserved':d['regions_val']['divisions'],'original_five_fold_assignment_available':False,'status':'configuration template only; not an exact archived manuscript configuration'})
 (DEST/'template_audit.json').write_text(json.dumps(receipt,indent=2)+'\n')
 print('Prepared four manuscript-setting templates. Original five-fold manifests remain required.')

if __name__=='__main__':main()
