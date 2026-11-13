from pathlib import Path
import os,sys,json,hashlib,subprocess
e=Path(os.environ['EVIDENCE_DIR']);root=Path(os.environ['GITHUB_WORKSPACE'])/'candidate'
mode=sys.argv[1]
sha=lambda b:hashlib.sha256(b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=root).decode().strip()
def save(name,obj):(e/name).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def inventory(p):return {f.relative_to(p).as_posix():sha(f.read_bytes()) for f in sorted(p.rglob('*')) if f.is_file()}
files={
 'tests/testthat/helper-issue199-portability.R':('502370eac752617963c987588fd6f570d9d4ce66','70772a1ead870aff71e8c505e6673915fd86fef6d7927005b069430360b65f48'),
 'tests/testthat/test-issue199-posterior-links.R':('1890b2ecdc414b936813d33c684e7cf32f2cccdf','73ab617b999e967d687719e1a041ced992fe680c5e5e5a3d96691d40005a8800'),
 'tests/testthat/test-shiny-example-parity.R':('7762365186d68eee350df5931179908259233260','722fc7f3aa394042e834104e88c993523c31ad0832c26d7725f85b8bf787bbae')}
result={'mode':mode,'verified':False}
try:
    result.update(commit=git('rev-parse','HEAD'),tree=git('rev-parse','HEAD^{tree}'),parent=git('show','-s','--format=%P'),status=git('status','--porcelain','--untracked-files=all'))
    save('source/identity-'+mode+'.json',result)
    assert result['commit']=='2e1120048a18d5413061e5cc6d4bf332f7fe42a7'
    assert result['tree']=='497a60ed3b9a64535a86bbaa6d1dcbf7fbf8aa58'
    assert result['parent']=='5bcb37d9f90b7f8d62949432fdeb867f951eb5e5'
    assert result['status']==''
    assert git('diff-tree','--no-commit-id','--name-only','-r','HEAD').splitlines()==list(files)[:2]
    for f,(blob,digest) in files.items():
        assert git('rev-parse','HEAD:'+f)==blob
        assert sha((root/f).read_bytes())==digest
    b1=list(files)[2]
    assert git('rev-parse','HEAD:'+b1)==git('rev-parse','HEAD^:'+b1)
    rows={}
    for item in subprocess.check_output(['git','ls-tree','-r','-z','HEAD'],cwd=root).split(b'\0'):
        if not item:continue
        meta,path=item.split(b'\t',1);f=path.decode();file_mode,kind,blob=meta.decode().split()
        content=(root/f).read_bytes()
        assert content==subprocess.check_output(['git','cat-file','blob',blob],cwd=root)
        rows[f]={'mode':file_mode,'blob':blob,'sha256':sha(content)}
    save('source/tracked-'+mode+'.json',rows)
    result['source_files']=len(rows)
    if mode!='before':
        assert rows==json.loads((e/'source/tracked-before.json').read_text())
        pkg=Path(os.environ['PACKAGE_LIB'])/'BreastSubtypeR'
        assert pkg.is_dir()
        installed=inventory(pkg);save('installed/files-'+mode+'.json',installed)
        tests={f:v['sha256'] for f,v in rows.items() if f.startswith('tests/') and f.endswith('.R')}
        for f,h in tests.items():assert installed[f]==h
        result.update(installed_files=len(installed),installed_R_test_files=len(tests))
        if mode=='after': assert installed==json.loads((e/'installed/files-installed.json').read_text())
    result['verified']=True
except Exception as exc:
    result['error']=repr(exc)
finally:
    save('source/identity-'+mode+'.json',result)
print(json.dumps(result,indent=2))
if not result['verified']:sys.exit(1)
