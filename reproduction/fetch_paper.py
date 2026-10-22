#!/usr/bin/env python3
"""Fetch the pinned, unmodified manuscript; reject a different PDF version."""
import hashlib,json,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def main():
 paper=json.loads((ROOT/'reproduction/paper_protocol/protocol.json').read_text())['paper']
 dest=ROOT/paper['local_pdf']
 if dest.exists(): data=dest.read_bytes()
 else:
  req=urllib.request.Request(paper['pdf_url'],headers={'User-Agent':'Mozilla/5.0'})
  with urllib.request.urlopen(req,timeout=90) as response: data=response.read()
  if not data.startswith(b'%PDF'): raise ValueError('Official URL did not return a PDF')
  if hashlib.sha256(data).hexdigest()!=paper['pdf_sha256']: raise ValueError('PDF differs from the pinned manuscript; review it before updating the protocol')
  dest.parent.mkdir(parents=True,exist_ok=True);dest.write_bytes(data)
 if hashlib.sha256(data).hexdigest()!=paper['pdf_sha256']: raise ValueError('Local manuscript checksum mismatch')
 print('Verified GHIST+ v1 manuscript:',dest)

if __name__=='__main__':main()
