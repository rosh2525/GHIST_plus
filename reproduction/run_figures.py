#!/usr/bin/env python3
"""Execute the authors' notebooks, preserve their tables, and export plots."""
from __future__ import annotations
import argparse, base64, datetime, hashlib, json, os, platform, subprocess, sys, time
from pathlib import Path

ROOT=Path(__file__).resolve().parents[1]
RESULTS=ROOT/'results'/'reproduction'
UPSTREAM='c5c895174c19903cc27d2661934b9b726bdb9306'

EXTRA={
2:"""
import shutil
from pathlib import Path
_persist = Path(os.environ['GHIST_REPRO_OUTPUT'])
for _src in OUT.iterdir():
    if _src.is_file(): shutil.copy2(_src, _persist / _src.name)
genes.to_csv(_persist / 'gene_summary.csv', index=False)
correlations.to_csv(_persist / 'pathway_correlations.csv', index=False)
pathways.to_csv(_persist / 'pathway_summary.csv', index=False)
""",
3:"""
import shutil
_persist = Path(os.environ['GHIST_REPRO_OUTPUT'])
for _src in WORK.iterdir():
    if _src.is_file(): shutil.copy2(_src, _persist / _src.name)
""",
4:"""
import shutil
_persist = Path(os.environ['GHIST_REPRO_OUTPUT'])
for _key, _value in list(globals().items()):
    if _key.startswith('fig4') and isinstance(_value, pd.DataFrame):
        _value.to_csv(_persist / (_key + '.csv'), index=False)
for _src in FIG4F_OUTPUT_DIR.iterdir():
    if _src.is_file(): shutil.copy2(_src, _persist / _src.name)
""",
5:"""
_persist = Path(os.environ['GHIST_REPRO_OUTPUT'])
paired.to_csv(_persist / 'paired_gene_pcc.csv', index=False)
coverage.to_csv(_persist / 'coverage_summary.csv', index=False)
_final = plot_figure5(paired)
_final.savefig(_persist / 'Figure5.png', dpi=300, facecolor='white')
_final.savefig(_persist / 'Figure5.pdf', facecolor='white')
plt.close(_final)
""",
}

def execute(number, core_only):
 import nbformat
 from nbclient import NotebookClient
 dest=RESULTS/f'Figure{number}'; dest.mkdir(parents=True,exist_ok=True)
 source=ROOT/f'Figure{number}.ipynb'; original=source.read_bytes()
 nb=nbformat.reads(original.decode(),as_version=4)
 for c in nb.cells:
  if c.cell_type=='code': c.outputs=[]; c.execution_count=None
 modifications=[]
 if number==3 and core_only:
  old='for panel, dataset in (("3h", "ATERA"), ("3i", "Breast5K")):'
  if old not in nb.cells[3].source: raise RuntimeError('Unexpected upstream Figure3 execution cell')
  nb.cells[3].source=nb.cells[3].source.replace(old,'for panel, dataset in ():')
  for line in ['    "3h": (plot_external(tables["3h"], "ATERA"), (4, 4)),\n','    "3i": (plot_external(tables["3i"], "Breast5K"), (4, 4)),\n']:
   if line not in nb.cells[4].source: raise RuntimeError('Unexpected upstream Figure3 plot cell')
   nb.cells[4].source=nb.cells[4].source.replace(line,'')
  modifications.append('Explicitly skip external panels 3h/3i; core panels and metric code are unchanged.')
 if number==2:
  # Figure2 does not import os. Added persistence cell requires it.
  nb.cells.append(nbformat.v4.new_code_cell('import os\n'+EXTRA[number]))
 else: nb.cells.append(nbformat.v4.new_code_cell(EXTRA[number]))
 modifications.append('Append output persistence cell; no model predictions, thresholds or metric functions changed.')
 os.environ['GHIST_REPRO_OUTPUT']=str(dest)
 print(f'Executing Figure{number}: {len(nb.cells)} cells',flush=True)
 started=time.monotonic()
 metadata={'figure':number,'upstream_commit':UPSTREAM,'source_notebook_sha256':hashlib.sha256(original).hexdigest(),'started_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'modifications':modifications,'model_training_performed':False,'new_model_inference_performed':False}
 client=NotebookClient(nb,timeout=7200,kernel_name='ghist-plus-repro',resources={'metadata':{'path':str(ROOT)}},allow_errors=False)
 def progress(cell,cell_index,**kwargs): print(f'Figure{number} cell {cell_index+1}/{len(nb.cells)}',flush=True)
 client.on_cell_start=progress
 try:
  client.execute()
  metadata['status']='passed'
 except Exception as exc:
  metadata['status']='failed'; metadata['error']=str(exc)
  raise
 finally:
  metadata['elapsed_seconds']=round(time.monotonic()-started,3)
  nbformat.write(nb,dest/f'Figure{number}.executed.ipynb')
  (dest/'execution.json').write_text(json.dumps(metadata,indent=2)+'\n')
  plot_count=0
  for ci,cell in enumerate(nb.cells):
   for oi,out in enumerate(cell.get('outputs',[])):
    png=out.get('data',{}).get('image/png')
    if png:
     plot_count+=1
     (dest/f'cell_{ci:02d}_output_{oi:02d}.png').write_bytes(base64.b64decode(''.join(png) if isinstance(png,list) else png))
  print(f'Figure{number}: {metadata["status"]}; {plot_count} displayed plots saved',flush=True)
 return metadata

def main():
 p=argparse.ArgumentParser(); p.add_argument('--figures',nargs='+',type=int,choices=[2,3,4,5],default=[2,4,5,3]); p.add_argument('--core-only',action='store_true'); args=p.parse_args()
 RESULTS.mkdir(parents=True,exist_ok=True)
 jupyter=RESULTS/'.jupyter'; kernel=jupyter/'kernels'/'ghist-plus-repro'; kernel.mkdir(parents=True,exist_ok=True)
 (kernel/'kernel.json').write_text(json.dumps({'argv':[sys.executable,'-m','ipykernel_launcher','-f','{connection_file}'],'display_name':'GHIST+ reproduction','language':'python'}))
 os.environ['JUPYTER_PATH']=str(jupyter)
 os.environ['JUPYTER_RUNTIME_DIR']=str(RESULTS/'.jupyter-runtime')
 os.environ['IPYTHONDIR']=str(RESULTS/'.ipython')
 os.environ['MPLCONFIGDIR']=str(RESULTS/'.matplotlib')
 os.environ['GHIST_BUNDLE_ROOT']=str(ROOT/'local_artifacts'/'bundle')
 os.environ['OMP_NUM_THREADS']='2'; os.environ['OPENBLAS_NUM_THREADS']='2'; os.environ['VECLIB_MAXIMUM_THREADS']='2'
 env={'python':sys.version,'platform':platform.platform(),'machine':platform.machine(),'executable':sys.executable,'upstream_commit':UPSTREAM,'scope':'released prediction evaluation and published-summary replotting; no new image inference or training'}
 (RESULTS/'environment.json').write_text(json.dumps(env,indent=2)+'\n')
 with (RESULTS/'pip_freeze.txt').open('w') as f: subprocess.run([sys.executable,'-m','pip','freeze'],stdout=f,check=True)
 manifests=[]
 for number in args.figures: manifests.append(execute(number,args.core_only))
 (RESULTS/'run_manifest.json').write_text(json.dumps(manifests,indent=2)+'\n')
 print('Selected notebook executions passed.',flush=True)

if __name__=='__main__': main()
