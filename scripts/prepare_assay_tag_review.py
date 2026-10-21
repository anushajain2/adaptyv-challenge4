"""Archive the documented platform-default suffix and prepare review inputs only."""
import json
import re
from datetime import datetime, timezone
from prepare_controls import ROOT
from receptor_preparation_v2 import sha

SUFFIX = 'GGGSGGGSRDHMVLHEYVNAAGITGGGSGGGSGGSAWSHPQFEKGGGSGGGSGGSAWSHPQFEK'
SOURCE_URL = 'https://docs.adaptyvbio.com/wiki/technologies/construct-tag-design'


def main():
    source = ROOT/'outputs/adaptyv-construct-tag-design.md'
    assert re.search(r'^'+SUFFIX+r'$',source.read_text(),re.M), 'Published suffix does not match'
    out = ROOT/'data/assay_tag_review'; out.mkdir(exist_ok=False)
    manifest = json.loads((ROOT/'data/expansion_his_confirmation/manifest.json').read_text())
    cases = []
    for c in manifest['cases']:
        bare = c['binder_sequence']
        cases.append({'id':c['id'],'kind':c['kind'],'parent_id':c.get('parent_id'),'bare_sequence':bare,
                      'tagged_sequence':bare+SUFFIX,'bare_length':len(bare),'tagged_length':len(bare+SUFFIX),
                      'suffix_start_1based':len(bare)+1,'suffix_end_1based':len(bare+SUFFIX)})
    record = {'source_url':SOURCE_URL,'retrieved_utc':datetime.now(timezone.utc).isoformat(),'source_sha256':sha(source),
        'suffix':SUFFIX,'suffix_length':len(SUFFIX),'suffix_histidine_positions_1based':[i+1 for i,x in enumerate(SUFFIX) if x=='H'],
        'evidence_status':'Exact publicly documented platform-default suffix; challenge page confirms C-terminal Twin-Strep capture but does not independently specify this full suffix.',
        'challenge_url':'https://proteinbase.com/competitions/anthropic-adaptyv-2026/challenges/tfr1',
        'facts':['Default C-terminal suffix includes flexible linkers, split-GFP11 and Twin-Strep.',
                 'The standard documentation says purification tags are not cleaved.',
                 'Designs should be submitted as bare sequences; these tagged strings are modeling inputs only.'],
        'prediction_status':'Not run; these files are review inputs, not submission sequences or evidence of tag compatibility.',
        'next_test':'Paired bare and tagged predictions with matched parents and controls, fixed seeds and unchanged binder-core metrics; inspect tag-receptor contacts separately. Do not count tag contacts as designed binding.',
        'source_manifest_sha256':sha(ROOT/'data/expansion_his_confirmation/manifest.json'),'cases':cases}
    (out/'manifest.json').write_text(json.dumps(record,indent=2))
    (out/'tagged_for_modeling_only.fasta').write_text(''.join(f">{c['id']} PLATFORM_DEFAULT_TAG_MODELING_ONLY\n{c['tagged_sequence']}\n" for c in cases))
    print(len(cases),'review constructs, suffix length',len(SUFFIX))

if __name__=='__main__':main()
