import datetime,hashlib,html,json,os,pathlib,subprocess
root=pathlib.Path('/Users/roshank/Documents/Codex/2026-10-02/we-x20')
repo=root/'GHIST_plus';out=root/'work/ghist_workflow_deck/evidence'
env=os.environ.copy();env['PATH']=str(root/'work/ghist_plus_reproduction/.venv/bin')+':'+env['PATH']
records=[]
for mode in ['data','checks','compare','status']:
 args=['python','reproduction/presentation_evidence.py',mode]
 r=subprocess.run(args,cwd=repo,env=env,capture_output=True,text=True)
 if r.returncode: raise RuntimeError(r.stdout+r.stderr)
 record={'mode':mode,'command':' '.join(args),'stdout':r.stdout,'stderr':r.stderr,'returncode':r.returncode,'captured_utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'script_sha256':hashlib.sha256((repo/'reproduction/presentation_evidence.py').read_bytes()).hexdigest()}
 records.append(record);(out/(mode+'.txt')).write_text('$ '+record['command']+'\n'+r.stdout)
 print(mode+'\n'+r.stdout,flush=True)
(out/'transcripts.json').write_text(json.dumps(records,indent=2)+'\n')
head='''<!doctype html><meta charset="utf-8"><title>GHIST+ command evidence</title><style>*{box-sizing:border-box}body{margin:0;background:#faf9f6;font-family:Arial,sans-serif}#capture{width:1160px;background:#102326;color:#f5f8f7;padding:20px 26px 22px}header{color:#a5c9c5;font-size:18px;margin-bottom:15px}pre{margin:0;font-family:Menlo,monospace;font-size:25px;line-height:1.42;white-space:pre-wrap;overflow-wrap:anywhere}.command{color:#8cdbca;margin-bottom:10px}footer{font-size:15px;color:#b9cecc;margin-top:14px}</style><main id="capture">'''
for record in records:
 body=head+'<header>GHIST+ / actual command transcript / 22 October 2026</header><pre class="command">$ '+html.escape(record['command'])+'</pre><pre>'+html.escape(record['stdout'].rstrip())+'</pre><footer>Browser capture of executed CLI output. The notebook results refer to the recorded 15 October run.</footer></main>'
 (out/(record['mode']+'.html')).write_text(body)
