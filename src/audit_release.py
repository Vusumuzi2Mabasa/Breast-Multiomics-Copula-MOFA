"""Verify the released benchmark, frozen sources and manifest without refitting."""
from pathlib import Path
import hashlib,json,sys
ROOT=Path(__file__).resolve().parents[1]


def main():
    lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
    for name,digest in lock['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    audit=json.loads((ROOT/'results/audit.json').read_text())
    assert audit['evaluation_datasets']==360 and audit['method_records']==2160
    assert audit['failed_fits']==0
    for name,digest in audit['cohort_hashes'].items():
        assert hashlib.sha256((ROOT/'results/replicates'/name).read_bytes()).hexdigest()==digest,name
    assert json.loads((ROOT/'results/validation/replay.json').read_text())['status']=='passed'
    assert json.loads((ROOT/'results/validation/R_metrics.json').read_text())['status']=='passed'
    doc=json.loads((ROOT/'results/document_structure.json').read_text())
    assert doc['author']=='Vusumuzi Mabasa' and doc['figures']==7 and doc['tables']==4 and doc['editable_equations']>=8
    manifest=ROOT/'checksums.sha256'
    if manifest.exists():
        entries=0
        for line in manifest.read_text().splitlines():
            digest,name=line.split('  ',1)
            assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
            entries+=1
        print('Manifest verified:',entries,'files')
    print('Release audit passed: 360 cohorts; 2160 records; validated sources, report and arithmetic.')


if __name__=='__main__':main()
