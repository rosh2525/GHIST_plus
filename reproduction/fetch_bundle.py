#!/usr/bin/env python3
"""Download only official GHIST+ figure inputs; pin and verify every file."""
from __future__ import annotations
import argparse, concurrent.futures, hashlib, json, shutil, time, urllib.request
from pathlib import Path

REPO = 'SydneyBioX/GHIST-Plus-bundle'
ROOT = Path(__file__).resolve().parents[1]
STATE = ROOT / 'reproduction'
BUNDLE = ROOT / 'local_artifacts' / 'bundle'
CORE_EXACT = {
 'evaluation_data/breast2/cell_gene_matrix_filtered.csv',
 'evaluation_data/breast2/cell_coords.csv',
 'evaluation_data/breast2/cell_type_filtered.csv',
 'evaluation_data/imputation/figure3_coordinates.csv.gz',
 'evaluation_data/imputation/figure3_ground_truth.csv.gz',
 'GHIST_plus/models/imputation/holdout_genes_by_slide.json',
 'GHIST_plus/predictions/breast_single/trainval_slide3_pred_expr_scaled.csv',
 'GHIST_plus/predictions/breast_multi/test_slide3_pred_expr_scaled.csv',
 'GHIST_plus/predictions/vq_off/vq_composition_off_trainval_slide3_pred_expr_scaled.csv',
 'GHIST_plus/predictions/vq_off/randomized_vq_composition_trainval_slide3_pred_expr_scaled.csv',
 'other_models/within_slide/GHIST_pred_expr.csv',
 'other_models/within_slide/SpatialEx_UNI2H_pred_expr.csv',
 'other_models/across_slide/GHIST_pred_expr.csv',
 'other_models/across_slide/SpatialEx_UNI2H_pred_expr.csv',
 'other_models/across_slide/Phoenix_pred_expr.csv',
 'README.md', 'THIRD_PARTY_NOTICES.md',
}
EXTERNAL_EXACT = {
 'GHIST_plus/predictions/breast_multi/atera/test_slide9101_pred_expr_scaled.npz',
 'GHIST_plus/predictions/breast_multi/breast5k/test_slide6_pred_expr_scaled.npz',
 'evaluation_data/atera/cell_coords.csv',
 'evaluation_data/atera/cell_gene_matrix_filtered.csv',
 'evaluation_data/breast5k/cell_coords_histology_filtered.csv',
 'evaluation_data/breast5k/cell_gene_matrix_filtered.csv',
}
CORE_PREFIXES = ['GHIST_plus/predictions/imputation/', 'other_models/imputation/', 'figure_data/figure5/']

def request_json(url):
 with urllib.request.urlopen(urllib.request.Request(url, headers={'User-Agent':'GHIST-local-reproduction/1.0'}), timeout=120) as r:
  return json.load(r)

def hashes(p):
 sha256=hashlib.sha256(); blob=hashlib.sha1(f'blob {p.stat().st_size}\0'.encode())
 with p.open('rb') as f:
  while chunk:=f.read(8*1024*1024): sha256.update(chunk); blob.update(chunk)
 return sha256.hexdigest(),blob.hexdigest()

def verify(p,item):
 if not p.is_file() or p.stat().st_size != item['size']: return None
 sha,blob=hashes(p)
 lfs=item.get('lfs') or {}
 expected=lfs.get('oid')
 if expected and sha!=expected: return None
 if not expected and blob!=item['oid']: return None
 return sha

def download(item,revision):
 relative=item['path']; dest=BUNDLE/relative; dest.parent.mkdir(parents=True,exist_ok=True)
 existing=verify(dest,item)
 if existing:
  print('verified cache',relative,flush=True)
  return {'path':relative,'bytes':item['size'],'sha256':existing,'source_verified':True}
 part=dest.with_name(dest.name+'.partial')
 started=time.monotonic()
 for attempt in range(1,4):
  offset=part.stat().st_size if part.exists() else 0
  if offset>item['size']: part.unlink(); offset=0
  if offset==item['size']:
   part.replace(dest)
   checksum=verify(dest,item)
   if checksum: return {'path':relative,'bytes':item['size'],'sha256':checksum,'source_verified':True}
   dest.unlink(); offset=0
  url=f'https://huggingface.co/datasets/{REPO}/resolve/{revision}/{relative}?download=true&repro_attempt={attempt}'
  headers={'User-Agent':'GHIST-local-reproduction/1.0'}
  if offset: headers['Range']=f'bytes={offset}-'
  print(f'download {relative}: {offset/1e6:.1f}/{item["size"]/1e6:.1f} MB',flush=True)
  try:
   with urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=120) as response:
    append=bool(offset and response.status==206)
    done=offset if append else 0; mark=time.monotonic()
    with part.open('ab' if append else 'wb') as f:
     while chunk:=response.read(4*1024*1024):
      f.write(chunk); done+=len(chunk)
      if time.monotonic()-mark>=20:
       print(f'progress {relative}: {done/1e6:.1f}/{item["size"]/1e6:.1f} MB',flush=True); mark=time.monotonic()
   part.replace(dest)
   checksum=verify(dest,item)
   if checksum is None: raise ValueError('size or upstream checksum mismatch')
   print(f'completed {relative} ({time.monotonic()-started:.1f}s)',flush=True)
   return {'path':relative,'bytes':item['size'],'sha256':checksum,'source_verified':True}
  except Exception as exc:
   print(f'retry {attempt}: {relative}: {type(exc).__name__}: {exc}',flush=True)
   if dest.exists(): dest.unlink()
   if attempt==3: raise
 raise RuntimeError(relative)

def main():
 parser=argparse.ArgumentParser(); parser.add_argument('--scope',choices=['core','all'],default='all'); parser.add_argument('--workers',type=int,default=3)
 args=parser.parse_args(); STATE.mkdir(parents=True,exist_ok=True)
 manifest_file=STATE/'data_manifest.json'
 if manifest_file.exists(): manifest=json.loads(manifest_file.read_text()); revision=manifest['revision']; all_items=manifest['upstream_files']
 else:
  revision=request_json(f'https://huggingface.co/api/datasets/{REPO}')['sha']
  all_items=request_json(f'https://huggingface.co/api/datasets/{REPO}/tree/{revision}?recursive=true&expand=false&limit=1000')
  manifest={'repository':REPO,'revision':revision,'upstream_files':all_items,'downloads':[]}; manifest_file.write_text(json.dumps(manifest,indent=2)+'\n')
 chosen=[i for i in all_items if i['type']=='file' and (i['path'] in CORE_EXACT or any(i['path'].startswith(p) for p in CORE_PREFIXES) or (args.scope=='all' and i['path'] in EXTERNAL_EXACT))]
 BUNDLE.mkdir(parents=True,exist_ok=True)
 remaining=sum(i['size'] for i in chosen if not (BUNDLE/i['path']).exists()); free=shutil.disk_usage(BUNDLE).free
 if remaining+5*1024**3>free: raise RuntimeError(f'Need {remaining/1e9:.2f} GB plus 5 GiB reserve; only {free/1e9:.2f} GB free')
 print(f'Pinned dataset {revision}; {len(chosen)} files, {sum(i["size"] for i in chosen)/1e9:.3f} GB',flush=True)
 previous={i['path']:i for i in manifest['downloads']}
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
  futures=[pool.submit(download,item,revision) for item in chosen]
  for future in concurrent.futures.as_completed(futures):
   entry=future.result(); previous[entry['path']]=entry
   manifest['downloads']=sorted(previous.values(),key=lambda i:i['path']); manifest['scope']=args.scope
   manifest_file.write_text(json.dumps(manifest,indent=2)+'\n')
 print('All selected inputs verified.',flush=True)

if __name__=='__main__': main()
