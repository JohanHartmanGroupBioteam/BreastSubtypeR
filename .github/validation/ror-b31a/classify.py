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
    conditions=rows(e/'tests/diagnostic-conditions.tsv')+rows(e/'tests/ror-conditions.tsv')
    classified=[]
    for x in conditions:
        matching=[r for r in failures if r['path'] in x['message'] and r['classification']=='established_B3_exact_historical_signature']
        known=((x['file'],x['test']) in [('test-expanded-examples.R','expanded current AUTO results reproduce the accepted complete objects'),('test-issue199-ror-bounds.R','ROR formula identity applies to fresh current-runtime results without changing RNG')] and
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
    assert len(diagnostic)==73 and len(set(x['file'] for x in diagnostic))==8
    assert all(summary[k] for k in ['primary_B2','secondary_B1','fresh_posteriors','mutations'])
    assert summary['unknown_numeric_leaves']==0 and summary['unknown_test_conditions']==0
    ror=rows(e/'tests/ror-blocks.tsv')
    routes=rows(e/'auto/ror-primary-routes.tsv')
    cells=rows(e/'auto/ror-primary-998.tsv')
    primary_checks=rows(e/'auto/ror-primary-checks.tsv')
    policy_visits=[r for dataset in ['OSLO2EMIT0obj','TCGABRCAobj'] for r in rows(e/'auto'/f'{dataset}-policy-visits.tsv')]
    assert len(policy_visits)==18 and all(x['visits']=='1' for x in policy_visits)
    summary['complete_example_policy_visits']=len(policy_visits)
    supplement=rows(e/'tests/supplement-controls.tsv')
    summary['primary_ROR']=len(routes)==10 and len(cells)==998 and sum(int(x['n']) for x in routes)==998 and sum(int(x['formula_identities']) for x in routes)==998 and all(x['within_ceiling']=='TRUE' and x['formula_identity']=='TRUE' for x in cells)
    summary['primary_ROR'] &= len(primary_checks)==2 and [int(x['eligible']) for x in primary_checks]==[502,496] and all(all(x[k]=='TRUE' for k in ['exact_routes','no_duplicate_missing_visits','policy_19_unchanged','RNG_preserved']) for x in primary_checks)
    summary['supplement']=len(supplement)==7 and all(x['detected']=='TRUE' for x in supplement)
    summary['ROR_blocks']=len(ror)
    summary['ROR_passes']=sum(int(x['passed']) for x in ror)
    summary['ROR_errors']=sum(x['error']=='TRUE' for x in ror)
    assert len(ror)==8 and summary['primary_ROR'] and summary['supplement']
    # Compare actual execution with the preserved exact-source Windows blocks.
    baseline={(x['file'],x['test']):x for x in rows(here/'windows-new-blocks.tsv')+rows(here/'windows-focused-blocks.tsv')}
    reconciled=[]
    for x in ror+diagnostic:
        b=baseline[(x['file'],x['test'])]
        attempted=int(x['passed'])+int(x['failed'])+(x['error']=='TRUE')
        interrupted=x['error']=='TRUE'
        not_reached=max(0,int(b['passed'])-attempted) if interrupted else 0
        reconciled.append(dict(file=x['file'],test=x['test'],windows_passed=b['passed'],native_passed=x['passed'],native_failed=x['failed'],native_error=x['error'],not_reached_after_error=not_reached,interrupted_comparison=int(interrupted)))
        if not interrupted: assert int(x['passed'])==int(b['passed']) and x['failed']=='0' and x['warning']=='0' and x['skipped']=='FALSE'
    write(e/'tests/block-reconciliation.tsv',reconciled,list(reconciled[0]))
    summary['not_reached_assertions']=sum(x['not_reached_after_error'] for x in reconciled)
    summary['interrupted_comparisons']=sum(x['interrupted_comparison'] for x in reconciled)
    assert len(failures)==22 and all(x['classification']=='established_B3_exact_historical_signature' for x in failures)
    assert sum('/cIHC.itr/' in x['path'] for x in failures)==18
    assert sum('/all.probs/' in x['path'] for x in failures)==4
    summary['acceptance']='PASS_WITH_CONDITIONS_PRIMARY_ROR_B1_B2_RECORDED_B3_FAILURES'

except Exception as exc:
    summary['error']=repr(exc)
(e/'acceptance.json').write_text(json.dumps(summary,indent=2)+'\n')
print(json.dumps(summary,indent=2))
if summary['acceptance']=='BLOCKED':sys.exit(1)
