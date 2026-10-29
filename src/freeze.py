"""Freeze design and numerical sources before evaluation; refuse overwrite."""
from pathlib import Path
import hashlib,json,datetime,importlib.metadata
ROOT=Path(__file__).resolve().parents[1]
target=ROOT/'protocol/source_lock.json'
if target.exists():raise SystemExit('Source lock already exists; create a versioned extension to change the design.')
files=['src/models.py','src/simulate.py','src/benchmark.py','src/dependence.py','src/run_all.py','protocol/design.json','tests/test_models.py']
lock=dict(frozen_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
          sources={f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in files},
          versions={p:importlib.metadata.version(p) for p in ['numpy','scipy','scikit-learn','mofapy2','pandas']})
target.write_text(json.dumps(lock,indent=2));print('Design and sources frozen.')
