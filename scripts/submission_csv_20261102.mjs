import fs from 'node:fs/promises';
import path from 'node:path';
import {createRequire} from 'node:module';
const require=createRequire(path.resolve('outputs/submission-authoring/package.json'));
const {Workbook}=await import(require.resolve('@oai/artifact-tool'));

const choice=process.argv[2];
if (!['original','revised'].includes(choice)) throw new Error('Choose original or revised explicitly');
const out=process.argv[3]||'reports/submission_2026-11-02';
const csv=await fs.readFile('reports/submission_draft/submission_draft.csv','utf8');
const workbook=await Workbook.fromCSV(csv,{sheetName:'Submission'});
const sheet=workbook.worksheets.getItem('Submission');
let rows=sheet.getRange('A2:C9').values;
const header=['name','sequence','molecule_class'];
if(JSON.stringify(sheet.getRange('A1:C1').values[0])!==JSON.stringify(header)) throw new Error('Unexpected source columns');
if(choice==='revised') {
  const m=JSON.parse(await fs.readFile('data/reserve_tag/manifest.json','utf8'));
  const c=m.cases.find(x=>x.id==='apical_short_008_A41H_mpnn03_bare');
  rows=[[c.base_id,c.core_sequence,'protein'],...rows.filter(x=>x[0]!=='expansion_apical_short_015')];
  sheet.getRange('A2:C9').values=rows;
}
if(rows.length!==8 || new Set(rows.map(x=>x[0])).size!==8 || new Set(rows.map(x=>x[1])).size!==8) throw new Error('Count/uniqueness failed');
for(const [name,sequence,kind] of rows) {
  if(!/^[A-Za-z0-9_]+$/.test(name)||!/^[ACDEFGHIKLMNPQRSTVWY]{10,250}$/.test(sequence)||kind!=='protein') throw new Error('Invalid submission row');
}
workbook.recalculate();
const check=await workbook.inspect({kind:'table',range:'Submission!A1:C9',include:'values,formulas',tableMaxRows:9,tableMaxCols:3,tableMaxCellChars:250,maxChars:5000});
await fs.writeFile('outputs/submission-authoring/csv_inspection.jsonl',check.ndjson);
const matrix=sheet.getRange('A1:C9').values;
// All fields are validated above to exclude CSV delimiters and quoting characters.
await fs.writeFile(`${out}/submission.csv`,matrix.map(r=>r.join(',')).join('\r\n')+'\r\n');
await fs.writeFile(`${out}/submission.fasta`,rows.map(([n,s])=>`>${n}\n${s}\n`).join(''));
await fs.writeFile('outputs/submission-authoring/selection.json',JSON.stringify({choice,rows:rows.map(([name,sequence,molecule_class],i)=>({rank:i+1,name,sequence,molecule_class}))},null,2));
console.log(`Validated eight bare sequences; ${choice} shortlist.`);
