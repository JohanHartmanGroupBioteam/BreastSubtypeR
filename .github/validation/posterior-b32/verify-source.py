from pathlib import Path
import os,sys,json,hashlib,subprocess
e=Path(os.environ['EVIDENCE_DIR']);root=Path(os.environ['GITHUB_WORKSPACE'])/'candidate'
mode=sys.argv[1]
sha=lambda b:hashlib.sha256(b).hexdigest()
def git(*args):return subprocess.check_output(['git',*args],cwd=root).decode().strip()
def save(name,obj):(e/name).write_text(json.dumps(obj,sort_keys=True,indent=2)+'\n')
def inventory(p):return {f.relative_to(p).as_posix():sha(f.read_bytes()) for f in sorted(p.rglob('*')) if f.is_file()}
files={'tests/testthat/helper-issue199-portability.R': ('473153f5ff1edce53a5dc910af7224db5352a458', 'b47acca636547b9ab4aa3624c4352c931bf596e57be6f436ce97fc6f23d12f4b'), 'tests/testthat/test-expanded-examples.R': ('0080a1a957cc59364cd13f2d43bfc4f3f825725c', '33ff00e85f1dc4d962c12f1faf1ca513a551fcc3e7be8efc7a4f64ffaac9d98d'), 'tests/testthat/test-issue199-posterior-links.R': ('a88dce5be217120b02091fe5c2bed87c5797cd17', 'c6e3e0f7e1cd7997a93fa3cda298c481ca06aa665b7f7b3db71d4a42714b320b'), 'tests/testthat/test-issue199-posterior-bounds.R': ('48592572ed4958bc0720c457d98b1a0969f45c8b', 'd048f2daf95669ceb72cb458104d14f858bd213abd8da666873ade044cb9be32'), 'tests/testthat/test-issue199-ror-bounds.R': ('20b542972b2a75a6da73375175ccbde29f0a14db', '515f322ae0eb60d2509b332b1287bc6ace244fb5f052cf797fcb84b2b7be4b12'), 'tests/testthat/test-shiny-example-parity.R': ('7762365186d68eee350df5931179908259233260', '722fc7f3aa394042e834104e88c993523c31ad0832c26d7725f85b8bf787bbae')}
result={'mode':mode,'verified':False}
try:
    result.update(commit=git('rev-parse','HEAD'),tree=git('rev-parse','HEAD^{tree}'),parent=git('show','-s','--format=%P'),status=git('status','--porcelain','--untracked-files=all'))
    save('source/identity-'+mode+'.json',result)
    assert result['commit']=='94762a40c55cf6dc4acd3d53e1882e1747721b7d'
    assert result['tree']=='957e2fa44d2daa7202732218c9b91670176fd858'
    assert result['parent']=='4dee0a3173599bcd8740e73b56a06305dd624b36'
    assert result['status']==''
    assert git('diff-tree','--no-commit-id','--name-only','-r','HEAD').splitlines()==sorted(list(files)[:4])
    if mode=='before':
        (e/'source/complete.patch').write_bytes(subprocess.check_output(['git','diff','--full-index','--binary','HEAD^','HEAD'],cwd=root))
        for f in list(files)[:4]:
            target=e/'source/reviewed'/f;target.parent.mkdir(parents=True,exist_ok=True)
            target.write_bytes((root/f).read_bytes())
    for f,(blob,digest) in files.items():
        assert git('rev-parse','HEAD:'+f)==blob
        assert sha((root/f).read_bytes())==digest
    for preserved in list(files)[4:]:
        assert git('rev-parse','HEAD:'+preserved)==git('rev-parse','HEAD^:'+preserved)
    b1=list(files)[5]
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
    fitted={}
    for row in [{'path': 'data/AIMSmodel.RData', 'sha256': '52b1190a09535d1df0ff09c06a4e0a180374da81649faaaf00194fd178b07b92', 'md5': '6c541625fe22667540ef7c11d9005290'}, {'path': 'data/sspbc.models.RData', 'sha256': '9cf0c03578c4ed7fc5a0b8b2c8a1eb683778df039e6b02d53ddc5d65d500634c', 'md5': '1c2ba126cc1ea8ecc15dd6dc4656285d'}, {'path': 'data/OSLO2EMIT0obj.RData', 'sha256': '1b7c1a6a2b1fb21f4968f185e6a64e90c0282ccdd7205aaf78a04cbdf75d6439', 'md5': '04612ddb44d51c8b9faef21c7de892a0'}, {'path': 'data/TCGABRCAobj.RData', 'sha256': 'b4ac2d8c898557453d69f5c7a3da81a1fbce4db4486bf0623adffcb1597f5d7d', 'md5': '67359af00cb83c2bc0eb2b990b82012a'}]:
        f=row['path']; raw=(root/f).read_bytes()
        assert sha(raw)==row['sha256'] and hashlib.md5(raw).hexdigest()==row['md5']
        if mode!='before': assert (pkg/f).read_bytes()==raw
        fitted[f]=row
    save('source/model-fixture-'+mode+'.json',fitted)
    result['verified']=True
except Exception as exc:
    result['error']=repr(exc)
finally:
    save('source/identity-'+mode+'.json',result)
print(json.dumps(result,indent=2))
if not result['verified']:sys.exit(1)
