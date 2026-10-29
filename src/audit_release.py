"""Verify the released benchmark, frozen sources and manifest without refitting."""
from pathlib import Path
import hashlib,json,sys,csv,math
ROOT=Path(__file__).resolve().parents[1]


def main():
    lock=json.loads((ROOT/'protocol/source_lock.json').read_text())
    for name,digest in lock['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    extension=json.loads((ROOT/'protocol/dependence_validation_lock.json').read_text())
    for name,digest in extension['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    correction=json.loads((ROOT/'protocol/metric_correction_lock.json').read_text())
    for name,digest in correction['sources'].items():
        assert hashlib.sha256((ROOT/name).read_bytes()).hexdigest()==digest,name
    audit=json.loads((ROOT/'results/audit.json').read_text())
    assert audit['evaluation_datasets']==360 and audit['method_records']==2160
    assert audit['failed_fits']==0
    evaluated=list(csv.DictReader((ROOT/'results/evaluation.csv').open()))
    assert len(evaluated)==2160
    for row in evaluated:
        assert row['status']=='ok'
        for metric in ['auc','ap','brier']:assert math.isfinite(float(row[metric]))
        if row['method']!='Clinical only':
            score=float(row['recovery']);old=float(row['original_rank_conditional_recovery'])
            assert 0<=score<=1+1e-10 and score<=old+1e-10
            assert int(float(row['estimated_rank']))<=3
    for name,digest in audit['cohort_hashes'].items():
        assert hashlib.sha256((ROOT/'results/replicates'/name).read_bytes()).hexdigest()==digest,name
    assert json.loads((ROOT/'results/validation/replay.json').read_text())['status']=='passed'
    assert json.loads((ROOT/'results/validation/R_metrics.json').read_text())['status']=='passed'
    doc=json.loads((ROOT/'results/document_structure.json').read_text())
    assert doc['author']=='Vusumuzi Mabasa' and doc['figures']==7 and doc['tables']==4 and doc['editable_equations']>=8
    assert doc['table_rows']==[7,7,4,9]
    assert 'nan%' not in (ROOT/'README.md').read_text().lower()
    qa=json.loads((ROOT/'results/report_qa.json').read_text())
    assert qa['status']=='passed' and qa['pages']==13
    assert hashlib.sha256((ROOT/'report/Copula_MOFA_Simulation_Report.docx').read_bytes()).hexdigest()==qa['word_sha256']
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
