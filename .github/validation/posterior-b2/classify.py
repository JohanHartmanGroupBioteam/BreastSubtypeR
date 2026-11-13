from pathlib import Path
import os,json,csv,sys
e=Path(os.environ['EVIDENCE_DIR']);here=Path(__file__).resolve().parent
def rows(p):
    with p.open(newline='',encoding='utf-8') as f:return list(csv.DictReader(f,delimiter='\t'))
def write(p,items,fields):
    with p.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=fields,delimiter='\t');w.writeheader();w.writerows(items)
summary={'acceptance':'BLOCKED','primary_B2':False,'secondary_B1':False,'fresh_posteriors':False,'mutations':False}
try:
    primary=rows(e/'tests/primary-blocks.tsv')
    summary['primary_B2']=len(primary)==9 and sum(int(x['passed']) for x in primary)==320 and all(x['failed']=='0' and x['error']=='FALSE' and x['warning']=='0' and x['skipped']=='FALSE' for x in primary)
    b1=rows(e/'b1/secondary-blocks.tsv')
    summary['secondary_B1']=len(b1)==14 and sum(int(x['passed']) for x in b1)==234 and sum(x['skipped']=='TRUE' for x in b1)==1 and all(x['failed']=='0' and x['error']=='FALSE' and x['warning']=='0' for x in b1)
    checks=rows(e/'auto/posterior-checks.tsv')
    summary['independent_posterior_checks']=len(checks)
    summary['fresh_posteriors']=bool(checks) and all(x['passed']=='TRUE' for x in checks) and not rows(e/'auto/structural-differences.tsv')
    decisions=[r for dataset in ['OSLO2EMIT0obj','TCGABRCAobj'] for r in rows(e/'auto'/f'{dataset}-decisions-blocks.tsv')]
    summary['fresh_posteriors'] &= len(decisions)==2 and all(x['failed']=='0' and x['error']=='FALSE' and x['warning']=='0' and x['skipped']=='FALSE' for x in decisions)
    mutations=rows(e/'mutations/evidence/mutation-summary.tsv')
    expected=['control','old_decision_gate','old_walker_gate','default_ties','numeric_equal_selection']
    summary['mutations']=len(mutations)==5 and [x['mutant'] for x in mutations]==expected and mutations[0]['passed']=='320' and all(x['warning']=='0' and x['skipped']=='0' for x in mutations) and all((int(x['failed'])+int(x['error'])==0 if x['mutant']=='control' else int(x['failed'])+int(x['error'])>0) for x in mutations)
    historic={(x['dataset'],x['path']):x for x in rows(here/'historical-leaves.tsv')}
    numeric=rows(e/'auto/numeric-leaves.tsv')
    failures=[]
    for x in numeric:
        if x['accepted']=='TRUE':continue
        old=historic.get((x['dataset'],x['path']))
        signature=['n','changed','policy','exceeded','actual_sha256','reference_sha256','masks','zeros','tiny']
        known=old is not None and old['accepted']=='FALSE' and all(x[k]==old[k] for k in signature)
        failures.append(dict(x,classification='established_B3_exact_historical_signature' if known else 'BLOCKED_unverified_provenance'))
    write(e/'auto/numerical-failure-inventory.tsv',failures,list(numeric[0])+['classification'])
    conditions=rows(e/'tests/diagnostic-conditions.tsv')
    classified=[]
    for x in conditions:
        matching=[r for r in failures if r['path'] in x['message'] and r['classification']=='established_B3_exact_historical_signature']
        known=(x['file']=='test-expanded-examples.R' and x['test']=='expanded current AUTO results reproduce the accepted complete objects' and
            ('unlisted value must be exact' in x['message'] or 'numeric bound exceeded' in x['message']) and bool(matching) and
            ('expectation_error' in x['class'] or 'expectation_failure' in x['class']))
        classified.append(dict(x,classification='established_B3' if known else 'BLOCKED_unverified_provenance',matched_paths=';'.join(sorted(set(r['dataset']+':'+r['path'] for r in matching)))))
    write(e/'tests/diagnostic-classified.tsv',classified,list(conditions[0])+['classification','matched_paths'] if conditions else ['file','test','index','class','message','classification','matched_paths'])
    diagnostic=rows(e/'tests/diagnostic-blocks.tsv')
    summary.update(diagnostic_blocks=len(diagnostic),diagnostic_passes=sum(int(x['passed']) for x in diagnostic),
        diagnostic_failed=sum(int(x['failed']) for x in diagnostic),diagnostic_errors=sum(x['error']=='TRUE' for x in diagnostic),
        diagnostic_warnings=sum(int(x['warning']) for x in diagnostic),diagnostic_skips=sum(x['skipped']=='TRUE' for x in diagnostic),
        B3_numeric_leaves=sum(x['classification']=='established_B3_exact_historical_signature' for x in failures),
        unknown_numeric_leaves=sum(x['classification'].startswith('BLOCKED') for x in failures),
        known_B3_test_conditions=sum(x['classification']=='established_B3' for x in classified),
        unknown_test_conditions=sum(x['classification'].startswith('BLOCKED') for x in classified))
    assert len(diagnostic)==33 and len(set(x['file'] for x in diagnostic))==4
    assert all(summary[k] for k in ['primary_B2','secondary_B1','fresh_posteriors','mutations'])
    assert summary['unknown_numeric_leaves']==0 and summary['unknown_test_conditions']==0
    summary['acceptance']='PASS_B2_B1_WITH_RECORDED_B3_FAILURES' if failures else 'PASS_B2_B1_NO_DIAGNOSTIC_FAILURES'
except Exception as exc:
    summary['error']=repr(exc)
(e/'acceptance.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
if summary['acceptance']=='BLOCKED':sys.exit(1)
