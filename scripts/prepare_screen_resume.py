"""Resume only predictions absent from the interrupted screen; never regenerate."""
import hashlib,json
from datetime import datetime,timezone
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]

def prepare():
    run=ROOT/'outputs/broad009_complete/broad-009-6wrw'
    checkpoint=json.loads((ROOT/'reports/broad_campaign/run_record/interruption_checkpoint.json').read_text())
    missing=set(checkpoint['missing_ids']);assert len(missing)==131 and not checkpoint['invalid_files']
    cases=[c for c in json.loads((run/'cases.json').read_text()) if c['id'] in missing]
    out=ROOT/'data/screen_resume';out.mkdir(exist_ok=False)
    for case in cases:
        for ext in ('.cif','.npz'):(out/(case['id']+ext)).write_bytes((run/'references'/(case['id']+ext)).read_bytes())
    m=json.loads((run/'inputs/manifest.json').read_text())
    m.update(cases=cases,source_run=run.name,max_prediction_seconds=1800,
        prepared_utc=datetime.now(timezone.utc).isoformat(),scope='Complete 131 missing screening predictions after generic Modal cancellation; no generation or reranking',
        random_state_note='Seed 307 restarted in a fresh process for missing cases only; random-number trajectory differs from an uninterrupted run. One sample per candidate remains the fixed screen.')
    m['input_sha256']={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file()}
    (out/'manifest.json').write_text(json.dumps(m,indent=2))
if __name__=='__main__':prepare()
