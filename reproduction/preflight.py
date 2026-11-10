#!/usr/bin/env python3
"""Report what is present/missing for fresh GHIST+ image inference or training."""
import argparse, importlib.util, json, platform
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
 p=argparse.ArgumentParser(); p.add_argument('--config',type=Path,default=ROOT/'configs/config_multi_breast.json'); p.add_argument('--data-root',type=Path,default=ROOT/'local_artifacts/training_data'); p.add_argument('--uni2-weights',type=Path); p.add_argument('--output',type=Path,default=ROOT/'results/reproduction/model_readiness.json'); a=p.parse_args()
 config=json.loads(a.config.read_text()); missing=[]; present=[]
 def visit(obj):
  if isinstance(obj,dict):
   for k,v in obj.items():
    if k.startswith('fp_') and isinstance(v,str):
     file=Path(v.replace('${GHIST_DATA_ROOT}',str(a.data_root)))
     record={'field':k,'path':str(file)}
     (present if file.exists() else missing).append(record)
    else: visit(v)
  elif isinstance(obj,list):
   for v in obj: visit(v)
 visit(config)
 torch_available=importlib.util.find_spec('torch') is not None
 cuda=False
 if torch_available:
  import torch
  cuda=bool(torch.cuda.is_available())
 local_encoder=bool(a.uni2_weights and a.uni2_weights.is_file())
 info={'config':str(a.config),'platform':platform.platform(),'machine':platform.machine(),'cuda_available':cuda,'pytorch_installed':torch_available,'uni2_weights_supplied':local_encoder,'present_input_fields':present,'missing_input_fields':missing,'ready_for_authors_cuda_model_run':bool(cuda and local_encoder and not missing),'note':'Figure analysis does not require model checkpoints or the UNI2-H encoder. Original image training/inference does. This check does not submit a gated-model access request or change the model.'}
 a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(json.dumps(info,indent=2)+'\n')
 print(f'CUDA: {cuda}; PyTorch installed: {torch_available}; UNI2-H supplied: {local_encoder}; missing input fields: {len(missing)}')
 print(f'Receipt: {a.output}')
 return 0 if info['ready_for_authors_cuda_model_run'] else 2

if __name__=='__main__': raise SystemExit(main())
