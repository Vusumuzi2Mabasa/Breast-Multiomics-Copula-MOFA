"""Check frozen sources, baseline preservation and all Study 3 outputs."""
import argparse,csv,hashlib,json,math
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];BASE=ROOT.parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--fresh',action='store_true');a=ap.parse_args()
    for name,prefix in [('dependence_lock.json',ROOT),('integration_lock_r2.json',BASE),('diablo_lock.json',BASE)]:
        lock=json.loads((ROOT/'protocol'/name).read_text())
        for path,digest in lock['sources'].items():assert sha(prefix/path)==digest,path
    preserved=0
    if not a.fresh:
        old=json.loads((BASE/'baseline/v1_preservation.json').read_text())
        for path,digest in old['sha256'].items():
            if path in old['allowed_root_updates']:continue
            assert sha(BASE/path)==digest,f'Original file changed: {path}'
            preserved+=1
    dep=json.loads((ROOT/'results/dependence_audit.json').read_text())
    assert dep['calibration_datasets']==5000 and dep['evaluation_datasets']==5800
    assert not dep['calibration_failures'] and not dep['evaluation_failures']
    assert dep['total_method_records']==34800
    assert sha(ROOT/'results/family_thresholds.json')==dep['threshold_sha256']
    cal_seeds=set();eval_seeds=set();chunks=0
    import numpy as np
    for phase,seen in [('calibration',cal_seeds),('evaluation',eval_seeds)]:
        for path in (ROOT/f'results/dependence_chunks/{phase}').glob('*.npz'):
            meta=json.loads(path.with_suffix('.json').read_text());assert sha(path)==meta['sha256']
            q=np.load(path);assert np.isfinite(q['statistics']).all();assert not meta['errors']
            seeds=set(q['seeds'].tolist());assert not seen&seeds;seen.update(seeds);chunks+=1
    assert not cal_seeds&eval_seeds
    ia=json.loads((ROOT/'results/integration_audit.json').read_text())
    assert ia['training_test_pairs']==96 and ia['method_records']==672 and ia['fit_failures']==0
    assert ia['diablo_fits']==96 and ia['diablo_failures']==0
    for file in ['replay.json','family_checks.json','R_conditional.json','integration_checks.json','R_integration.json']:
        assert json.loads((ROOT/'results/validation'/file).read_text())['status']=='passed',file
    assert len(list((ROOT/'figures').glob('*.png')))==7
    doc=json.loads((ROOT/'results/document_structure.json').read_text())
    assert doc['author']=='Vusumuzi Mabasa' and doc['figures']==7 and doc['tables']==6 and doc['editable_equations']>=10
    assert doc['table_rows']==[5,7,6,5,5,5]
    text=(BASE/'README.md').read_text();assert 'nan%' not in text.lower() and 'None%' not in text
    if not a.fresh:
        qa=json.loads((ROOT/'results/report_qa.json').read_text());assert qa['status']=='passed'
        assert sha(ROOT/'report/Study3_Coordination_Report_v2.docx')==qa['word_sha256']
    result=dict(status='passed',preserved_v1_files=preserved,dependence_chunks=chunks,calibration_datasets=len(cal_seeds),evaluation_datasets=len(eval_seeds),integration_replicates=96,diablo_replicates=96,word_pages=None if a.fresh else qa['pages'])
    if not a.fresh:(ROOT/'results/release_audit.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':main()
