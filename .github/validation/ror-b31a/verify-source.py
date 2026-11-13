from pathlib import Path
import os,sys,json,hashlib,subprocess
e=Path(os.environ['EVIDENCE_DIR']);root=Path(os.environ['GITHUB_WORKSPACE'])/'candidate'
mode=sys.argv[1]
sha=lambda b:hashlib.sha256(b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=root).decode().strip()
def save(name,obj):(e/name).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def inventory(p):return {f.relative_to(p).as_posix():sha(f.read_bytes()) for f in sorted(p.rglob('*')) if f.is_file()}
files={'tests/testthat/helper-issue199-portability.R': ('e54edb1a903c8ff6f3acbaae8e5efccb76238d6b', 'd6ac109502c57c68eca7354f0ed70da705fed4594755830366773f73ee33109f'), 'tests/testthat/test-expanded-examples.R': ('3f23abf535ded1829a6709ddb20a9e3bb43a61c7', 'f46a0e57f01cc63c10bf3d7d8e7f6d77ec9095d335d338dfe33c75bd77982ae2'), 'tests/testthat/test-issue199-ror-bounds.R': ('20b542972b2a75a6da73375175ccbde29f0a14db', '515f322ae0eb60d2509b332b1287bc6ace244fb5f052cf797fcb84b2b7be4b12'), 'tests/testthat/test-issue199-posterior-links.R': ('1890b2ecdc414b936813d33c684e7cf32f2cccdf', '73ab617b999e967d687719e1a041ced992fe680c5e5e5a3d96691d40005a8800'), 'tests/testthat/test-shiny-example-parity.R': ('7762365186d68eee350df5931179908259233260', '722fc7f3aa394042e834104e88c993523c31ad0832c26d7725f85b8bf787bbae')}
result={'mode':mode,'verified':False}
try:
    result.update(commit=git('rev-parse','HEAD'),tree=git('rev-parse','HEAD^{tree}'),parent=git('show','-s','--format=%P'),status=git('status','--porcelain','--untracked-files=all'))
    save('source/identity-'+mode+'.json',result)
    assert result['commit']=='7946ba5b8520584dec8de32fd83662951d1422a8'
    assert result['tree']=='254ec1ca636a561bd13b4909dddf6eda28fecfcf'
    assert result['parent']=='ce0414d1b35800985c1fa4b50957db1285eed848'
    assert result['status']==''
    assert git('diff-tree','--no-commit-id','--name-only','-r','HEAD').splitlines()==list(files)[:3]
    if mode=='before':
        (e/'source/complete.patch').write_bytes(subprocess.check_output(['git','diff','--full-index','--binary','HEAD^','HEAD'],cwd=root))
        for f in list(files)[:3]:
            target=e/'source/reviewed'/f;target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes((root/f).read_bytes())
    for f,(blob,digest) in files.items():
        assert git('rev-parse','HEAD:'+f)==blob
        assert sha((root/f).read_bytes())==digest
    for preserved in list(files)[3:]:
        assert git('rev-parse','HEAD:'+preserved)==git('rev-parse','HEAD^:'+preserved)
    b1=list(files)[4]
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
        for f,v in rows.items():
            if f.startswith('data/') or f.startswith('tests/testthat/fixtures/'):
                assert installed[f]==v['sha256'], f
        result.update(installed_files=len(installed),installed_R_test_files=len(tests))
        if mode=='after': assert installed==json.loads((e/'installed/files-installed.json').read_text())
    result['verified']=True
except Exception as exc:
    result['error']=repr(exc)
finally:
    save('source/identity-'+mode+'.json',result)
print(json.dumps(result,indent=2))
if not result['verified']:sys.exit(1)
