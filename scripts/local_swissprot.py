"""Local, checksum-verified fallback for the stalled NCBI Swiss-Prot job.

Download the official files named in prepare() into outputs/local_swissprot,
then run prepare and run. No network, paid compute, or system installation.
"""
import argparse
import gzip
import hashlib
import json
import shutil
import subprocess
import tarfile
import xml.etree.ElementTree as ET

from novelty_search import ROOT, OUT, now, save, sha

WORK = ROOT/'outputs/local_swissprot'
DEST = OUT/'local_swissprot'
BASE = 'https://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/'
BIN_URL = 'https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.17.0/ncbi-blast-2.17.0+-aarch64-macosx.tar.gz'


def digest(path, algorithm):
    with path.open('rb') as f:
        return hashlib.file_digest(f, algorithm).hexdigest()


def prepare():
    if (DEST/'plan.json').exists():
        raise SystemExit('Local fallback plan already frozen.')
    DEST.mkdir(parents=True, exist_ok=True)
    metadata=ET.parse(WORK/'RELEASE.metalink').getroot()
    ns={'m':'http://www.metalinker.org/'}
    f=next(f for f in metadata.findall('m:files/m:file',ns) if f.get('name')=='uniprot_sprot.fasta.gz')
    original=json.loads((OUT/'blast/request.json').read_text())
    (DEST/'queries.fasta').write_text(original['QUERY'])
    shutil.copyfile(WORK/'RELEASE.metalink',DEST/'RELEASE.metalink')
    shutil.copyfile(WORK/'blast.tar.gz.md5',DEST/'blast.tar.gz.md5')
    save(DEST/'plan.json',{
        'prepared_utc':now(), 'reason':'NCBI RID DHJRT8TH014 remains WAITING after more than one hour; local fallback avoids dependence on the remote queue.',
        'original_plan_sha256':sha(OUT/'plan.json'), 'query_sha256':sha(DEST/'queries.fasta'),
        'database_url':BASE+f.get('name'), 'metadata_url':BASE+'RELEASE.metalink',
        'release':metadata.findtext('m:version',namespaces=ns),
        'database_expected_bytes':int(f.findtext('m:size',namespaces=ns)),
        'database_expected_md5':f.findtext('m:verification/m:hash',namespaces=ns),
        'binary_url':BIN_URL, 'binary_md5_url':BIN_URL+'.md5',
        'binary_expected_md5':(WORK/'blast.tar.gz.md5').read_text().split()[0],
        'parameters':{'task':'blastp','evalue':10,'seg':'yes','soft_masking':'true',
                      'comp_based_stats':2,'matrix':'BLOSUM62','word_size':3,'max_target_seqs':100,
                      'num_threads':4,'outfmt':15},
        'scope':'Official UniProtKB/Swiss-Prot canonical FASTA, excluding the separate isoform supplement. This is a new explicitly versioned local search, not the completed result of the NCBI RID and not necessarily an identical database snapshot.',
        'interpretation':'Retain frozen identity/coverage/E-value heuristics per HSP. Human P02786 must be retrieved at >=99% identity and >=99% query coverage with E<=1e-10. Report masking and short-sequence limitations. No absence-of-homology or novel-fold claim.',
        'new_modal_compute_usd':0})
    print('Frozen local fallback plan for release',metadata.findtext('m:version',namespaces=ns))


def run():
    plan=json.loads((DEST/'plan.json').read_text())
    if (DEST/'results.json').exists():raise SystemExit('Results already exist; refusing overwrite.')
    assert sha(OUT/'plan.json')==plan['original_plan_sha256']
    assert sha(DEST/'queries.fasta')==plan['query_sha256']
    db=WORK/'uniprot_sprot.fasta.gz';binary=WORK/'blast.tar.gz'
    assert db.stat().st_size==plan['database_expected_bytes']
    assert digest(db,'md5')==plan['database_expected_md5']
    assert digest(binary,'md5')==plan['binary_expected_md5']
    with tarfile.open(binary) as tar:tar.extractall(WORK,filter='data')
    bindir=WORK/'ncbi-blast-2.17.0+'/ 'bin'
    blast=bindir/'blastp';make=bindir/'makeblastdb'
    version=subprocess.check_output([str(blast),'-version'],text=True)
    fasta=WORK/'swissprot.fasta'
    with gzip.open(db,'rb') as src, fasta.open('wb') as dst:shutil.copyfileobj(src,dst)
    records=0;letters=0
    with fasta.open() as f:
        for line in f:
            if line.startswith('>'):records+=1
            else:letters+=len(line.strip())
    dbpath=WORK/'swissprot'
    makecmd=[str(make),'-in',str(fasta),'-dbtype','prot','-parse_seqids','-out',str(dbpath),
             '-title',f"UniProtKB Swiss-Prot canonical {plan['release']}"]
    cmd=[str(blast),'-query',str(DEST/'queries.fasta'),'-db',str(dbpath),'-out',str(DEST/'results.tmp.json')]
    for k,v in plan['parameters'].items():cmd.extend(['-'+k,str(v)])
    record={'started_utc':now(),'state':'RUNNING','version':version,
            'database_sequences':records,'database_residues':letters,
            'database_archive_sha256':sha(db),'database_fasta_sha256':sha(fasta),
            'binary_archive_sha256':sha(binary),'blastp_sha256':sha(blast),
            'makeblastdb_command':makecmd,'blastp_command':cmd,'plan_sha256':sha(DEST/'plan.json')}
    save(DEST/'run.json',record)
    for label,args in [('makeblastdb',makecmd),('blastp',cmd)]:
        print('Running',label,flush=True)
        with (DEST/f'{label}.stdout.txt').open('w') as out,(DEST/f'{label}.stderr.txt').open('w') as err:
            result=subprocess.run(args,stdout=out,stderr=err)
        if result.returncode:
            record.update(state='FAILED',failed_step=label,returncode=result.returncode,finished_utc=now())
            save(DEST/'run.json',record);raise SystemExit(f'{label} failed: {result.returncode}')
    result=json.loads((DEST/'results.tmp.json').read_text())
    assert len(result['BlastOutput2'])==7
    (DEST/'results.tmp.json').rename(DEST/'results.json')
    record.update(state='COMPLETE',finished_utc=now(),result_sha256=sha(DEST/'results.json'))
    save(DEST/'run.json',record)
    print('Local search complete:',records,'sequences;',letters,'residues.')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('action',choices=['prepare','run'])
    globals()[p.parse_args().action]()
