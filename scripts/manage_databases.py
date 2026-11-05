#!/usr/bin/env python3
"""Managed metagenomic database utility v0.3.0.

Commands: list, info, status, verify, install.
Uses only the Python standard library.
"""
from __future__ import annotations
import argparse, hashlib, json, os, shutil, sys, tarfile, tempfile, time
import urllib.error, urllib.request
from pathlib import Path
from typing import Any

VERSION="0.3.0"; SCHEMA="1"; MANIFEST="installation.json"
OK=0; INSTALL_ERROR=1; USAGE_ERROR=2; UNKNOWN_DB=3; NO_SPACE=4; CONFLICT=5

class ManagedError(Exception): pass
class NoSpaceError(ManagedError): pass
class ConflictError(ManagedError): pass

def eprint(*x: object)->None: print(*x,file=sys.stderr)
def project_root()->Path: return Path(__file__).resolve().parent.parent

def abs_path(v:str, base:Path|None=None)->Path:
    p=Path(v).expanduser()
    return ((base or Path.cwd())/p if not p.is_absolute() else p).resolve(strict=False)

def load_catalog(path:Path)->dict[str,Any]:
    if path.is_dir(): raise ManagedError(f"--catalog expects a JSON file, not a directory: {path}")
    try: data=json.loads(path.read_text(encoding="utf-8"))
    except FileNotFoundError as ex: raise ManagedError(f"Catalog not found: {path}") from ex
    except (OSError,json.JSONDecodeError) as ex: raise ManagedError(f"Cannot read catalog '{path}': {ex}") from ex
    if not isinstance(data,dict) or data.get("schema_version")!=SCHEMA: raise ManagedError("Catalog must be an object with schema_version='1'.")
    dbs=data.get("databases")
    if not isinstance(dbs,dict) or not dbs: raise ManagedError("Catalog 'databases' must be a non-empty object.")
    default=data.get("default_database")
    if default not in dbs: raise ManagedError("default_database is absent from databases.")
    for dbid,e in dbs.items(): validate_entry(dbid,e)
    return data

def validate_entry(dbid:str,e:Any)->None:
    if not isinstance(e,dict): raise ManagedError(f"Database '{dbid}' must be an object.")
    req=("family","display_name","description","release","provider","archive_url","checksum_type","archive_checksum","archive_filename","archive_type","install_subdir","required_files")
    for k in req:
        if k not in e: raise ManagedError(f"Database '{dbid}' lacks '{k}'.")
    if e["checksum_type"] not in {"md5","sha256"}: raise ManagedError(f"Database '{dbid}' has unsupported checksum_type.")
    n=32 if e["checksum_type"]=="md5" else 64; h=e["archive_checksum"]
    if not isinstance(h,str) or len(h)!=n or any(c not in "0123456789abcdef" for c in h): raise ManagedError(f"Database '{dbid}' has an invalid checksum.")
    if not isinstance(e["required_files"],list) or not e["required_files"]: raise ManagedError(f"Database '{dbid}' required_files is invalid.")
    sub=Path(e["install_subdir"])
    if sub.is_absolute() or ".." in sub.parts: raise ManagedError(f"Database '{dbid}' install_subdir is unsafe.")
    b=e.get("bracken")
    if b:
        if not isinstance(b.get("read_lengths"),list) or not isinstance(b.get("required_file_template"),str): raise ManagedError(f"Database '{dbid}' Bracken metadata is invalid.")

def get_entry(cat:dict[str,Any],dbid:str|None)->tuple[str,dict[str,Any]]:
    dbid=dbid or cat["default_database"]
    if dbid not in cat["databases"]: raise KeyError(dbid)
    return dbid,cat["databases"][dbid]

def install_path(root:Path,e:dict[str,Any])->Path:
    p=(root/e["install_subdir"]).resolve(strict=False); rr=root.resolve(strict=False)
    try:p.relative_to(rr)
    except ValueError as ex: raise ManagedError(f"Unsafe installation path: {p}") from ex
    return p

def bracken_files(e:dict[str,Any])->list[str]:
    b=e.get("bracken") or {}
    if not b.get("included"): return []
    return [b["required_file_template"].format(read_length=x) for x in b["read_lengths"]]

def inspect(root:Path,dbid:str,e:dict[str,Any])->dict[str,Any]:
    p=install_path(root,e); m=p/MANIFEST
    r={"database_id":dbid,"release":e["release"],"path":str(p),"status":"not-installed","missing_files":[],"problems":[]}
    if not p.exists(): return r
    if not p.is_dir() or p.is_symlink(): r["status"]="conflict";r["problems"].append("Installation path is not a regular directory.");return r
    if not m.is_file() or m.is_symlink(): r["status"]="unmanaged";r["problems"].append("installation.json is missing or unsafe.");return r
    try: md=json.loads(m.read_text(encoding="utf-8"))
    except Exception as ex:r["status"]="invalid";r["problems"].append(f"Invalid installation.json: {ex}");return r
    for k,v in {"catalog_id":dbid,"release":e["release"],"status":"installed"}.items():
        if md.get(k)!=v:r["problems"].append(f"installation.json {k}={md.get(k)!r}; expected {v!r}.")
    for f in e["required_files"]+bracken_files(e):
        q=p/f
        if not q.is_file() or q.is_symlink() or q.stat().st_size<=0:r["missing_files"].append(f)
    r["status"]="installed" if not r["problems"] and not r["missing_files"] else "invalid";r["manifest"]=md
    return r

def verify(root:Path,dbid:str,e:dict[str,Any])->dict[str,Any]:
    r=inspect(root,dbid,e);r["verified"]=False;r["checked_files"]=[]
    if r["status"]!="installed":return r
    p=Path(r["path"]); md=r["manifest"]
    for k,v in {"checksum_type":e["checksum_type"],"archive_checksum":e["archive_checksum"],"archive_filename":e["archive_filename"]}.items():
        if k in md and md[k]!=v:r["problems"].append(f"installation.json {k} does not match catalog.")
    for f in e["required_files"]+bracken_files(e):
        q=p/f; good=q.is_file() and not q.is_symlink() and q.stat().st_size>0
        r["checked_files"].append({"path":f,"valid":good,"size_bytes":q.stat().st_size if good else None})
        if not good and f not in r["missing_files"]:r["missing_files"].append(f)
    r["verified"]=not r["problems"] and not r["missing_files"];r["status"]="verified" if r["verified"] else "invalid"
    return r

def hash_file(p:Path,algo:str)->str:
    d=hashlib.new(algo)
    with p.open("rb") as h:
        for chunk in iter(lambda:h.read(8*1024*1024),b""):d.update(chunk)
    return d.hexdigest()

def head_size(url:str)->int|None:
    try:
        req=urllib.request.Request(url,method="HEAD",headers={"User-Agent":f"manage-databases/{VERSION}"})
        with urllib.request.urlopen(req,timeout=30) as r:
            x=r.headers.get("Content-Length");return int(x) if x and x.isdigit() else None
    except Exception:return None

def ensure_space(root:Path,e:dict[str,Any],remote:int|None,part:Path)->None:
    installed=e.get("installed_size_bytes")
    if installed is None and e.get("reported_index_size_gb") is not None:installed=int(float(e["reported_index_size_gb"])*1_000_000_000)
    remain=max(0,remote-(part.stat().st_size if part.exists() else 0)) if remote is not None else None
    if remain is None or installed is None:eprint("WARNING: exact required space is unknown.");return
    required=remain+int(installed*1.10); parent=root
    while not parent.exists() and parent!=parent.parent:parent=parent.parent
    free=shutil.disk_usage(parent).free
    print(f"Disk space: approximately {required} bytes required; {free} available.")
    if free<required:raise NoSpaceError(f"Insufficient disk space: required {required}, available {free} bytes.")

def download(url:str,part:Path,expected:int|None)->None:
    part.parent.mkdir(parents=True,exist_ok=True);offset=part.stat().st_size if part.exists() else 0
    headers={"User-Agent":f"manage-databases/{VERSION}"}
    if offset:headers["Range"]=f"bytes={offset}-"
    try:r=urllib.request.urlopen(urllib.request.Request(url,headers=headers),timeout=60)
    except urllib.error.URLError as ex:raise ManagedError(f"Download failed: {ex}") from ex
    with r:
        append=offset>0 and getattr(r,"status",r.getcode())==206
        if offset and not append:eprint("WARNING: server ignored Range; restarting download.");offset=0
        with part.open("ab" if append else "wb") as out:
            total=offset;last=time.monotonic()
            while True:
                c=r.read(8*1024*1024)
                if not c:break
                out.write(c);total+=len(c)
                if time.monotonic()-last>=5:print(f"Downloaded {total}"+(f"/{expected} bytes" if expected else " bytes"));last=time.monotonic()
            out.flush();os.fsync(out.fileno())
    if expected is not None and part.stat().st_size!=expected:raise ManagedError(f"Downloaded size mismatch: {part.stat().st_size} != {expected}.")

def safe_extract(archive:Path,dest:Path)->None:
    dest.mkdir(mode=0o700,parents=False,exist_ok=False);root=dest.resolve()
    with tarfile.open(archive,"r:gz") as t:
        members=t.getmembers()
        if not members:raise ManagedError("Archive is empty.")
        for m in members:
            n=Path(m.name)
            if n.is_absolute() or ".." in n.parts or m.issym() or m.islnk() or m.isdev() or m.isfifo():raise ManagedError(f"Unsafe archive member: {m.name}")
            target=(root/n).resolve(strict=False)
            try:target.relative_to(root)
            except ValueError as ex:raise ManagedError(f"Archive path escapes staging: {m.name}") from ex
            if m.isdir():target.mkdir(parents=True,exist_ok=True)
            elif m.isfile():
                target.parent.mkdir(parents=True,exist_ok=True);src=t.extractfile(m)
                if src is None:raise ManagedError(f"Cannot read archive member: {m.name}")
                with src,target.open("xb") as out:shutil.copyfileobj(src,out,8*1024*1024)
            else:raise ManagedError(f"Unsupported archive member: {m.name}")

def locate_root(extracted:Path,e:dict[str,Any])->Path:
    names=e["required_files"]+bracken_files(e); candidates=[extracted]+[x for x in extracted.iterdir() if x.is_dir() and not x.is_symlink()]
    matches=[x for x in candidates if all((x/n).is_file() for n in names)]
    if len(matches)!=1:raise ManagedError("Could not identify one extracted database root with all required files.")
    return matches[0]

def atomic_json(p:Path,v:dict[str,Any])->None:
    p.parent.mkdir(parents=True,exist_ok=True); fd,name=tempfile.mkstemp(prefix=".manifest.",dir=p.parent)
    q=Path(name)
    try:
        with os.fdopen(fd,"w",encoding="utf-8") as h:json.dump(v,h,indent=2,sort_keys=True);h.write("\n");h.flush();os.fsync(h.fileno())
        os.replace(q,p)
    finally:q.unlink(missing_ok=True)

def install(root:Path,dbid:str,e:dict[str,Any])->dict[str,Any]:
    state=inspect(root,dbid,e)
    if state["status"]=="installed":return {"action":"already-installed","installation":verify(root,dbid,e)}
    if state["status"]!="not-installed":raise ConflictError(f"Refusing to overwrite {state['status']} installation: {state['path']}")
    downloads=root/"downloads";stages=root/".staging";target=install_path(root,e);archive=downloads/e["archive_filename"];part=downloads/(e["archive_filename"]+".part");stage=stages/f"{dbid}.{e['release']}.{os.getpid()}"
    if target.exists() or stage.exists():raise ConflictError("Installation target or staging directory already exists.")
    remote=e.get("archive_size_bytes") or head_size(e["archive_url"]);ensure_space(root,e,remote,part)
    downloads.mkdir(parents=True,mode=0o700,exist_ok=True);stages.mkdir(parents=True,mode=0o700,exist_ok=True);target.parent.mkdir(parents=True,mode=0o700,exist_ok=True)
    if archive.exists():
        if hash_file(archive,e["checksum_type"])!=e["archive_checksum"]:raise ConflictError(f"Existing archive checksum mismatch: {archive}")
    else:
        download(e["archive_url"],part,remote)
        actual=hash_file(part,e["checksum_type"])
        if actual!=e["archive_checksum"]:raise ManagedError(f"Checksum mismatch: {actual}; expected {e['archive_checksum']}. Partial retained: {part}")
        os.replace(part,archive)
    stage.mkdir(mode=0o700);extracted=stage/"extracted"
    try:
        safe_extract(archive,extracted);dbroot=locate_root(extracted,e)
        atomic_json(dbroot/MANIFEST,{"schema_version":"1","catalog_id":dbid,"family":e["family"],"release":e["release"],"status":"installed","installed_at":time.strftime("%Y-%m-%dT%H:%M:%SZ",time.gmtime()),"source_url":e["archive_url"],"archive_filename":e["archive_filename"],"checksum_type":e["checksum_type"],"archive_checksum":e["archive_checksum"],"required_files_verified":True})
        os.replace(dbroot,target)
    except Exception:shutil.rmtree(stage,ignore_errors=True);raise
    shutil.rmtree(stage,ignore_errors=True);result=verify(root,dbid,e)
    if not result["verified"]:raise ManagedError("Published installation failed verification.")
    return {"action":"installed","installation":result,"archive":str(archive)}

def human(n:Any)->str:
    if n is None:return "unknown"
    x=float(n)
    for u in ("B","KiB","MiB","GiB","TiB"):
        if x<1024 or u=="TiB":return f"{x:.1f} {u}"
        x/=1024
    return "unknown"
def emit(v:Any)->None:print(json.dumps(v,indent=2,sort_keys=True,ensure_ascii=False))

def context(a:argparse.Namespace)->tuple[dict[str,Any],Path]:
    pr=project_root();return load_catalog(abs_path(a.catalog,pr)),abs_path(a.root,pr)

def cmd_list(cat:dict[str,Any],root:Path,a:argparse.Namespace)->int:
    rows=[]
    for i in sorted(cat["databases"]):
        e=cat["databases"][i];st=inspect(root,i,e);rows.append({"id":i,"default":i==cat["default_database"],"family":e["family"],"release":e["release"],"status":st["status"]})
    if a.json:emit(rows)
    else:
        print("ID\tDEFAULT\tFAMILY\tRELEASE\tSTATUS")
        for r in rows:print(f"{r['id']}\t{'yes' if r['default'] else 'no'}\t{r['family']}\t{r['release']}\t{r['status']}")
    return OK

def cmd_info(cat:dict[str,Any],root:Path,a:argparse.Namespace)->int:
    i,e=get_entry(cat,a.database_id);r={"id":i,**e,"installation":inspect(root,i,e)}
    if a.json:emit(r)
    else:
        print(f"ID: {i}\nName: {e['display_name']}\nRelease: {e['release']}\nProvider: {e['provider']}\nArchive: {e['archive_filename']}\nInstallation: {r['installation']['status']}\nPath: {r['installation']['path']}")
    return OK

def cmd_status(cat:dict[str,Any],root:Path,a:argparse.Namespace)->int:
    ids=[a.database_id] if a.database_id else sorted(cat["databases"]);states=[]
    for x in ids:i,e=get_entry(cat,x);states.append(inspect(root,i,e))
    if a.json:emit(states)
    else:
        for r in states:print(f"{r['database_id']}@{r['release']}: {r['status']} ({r['path']})")
    return CONFLICT if any(r["status"] in {"invalid","conflict"} for r in states) else OK

def cmd_verify(cat:dict[str,Any],root:Path,a:argparse.Namespace)->int:
    i,e=get_entry(cat,a.database_id);r=verify(root,i,e)
    if a.json:emit(r)
    else:print(f"Database: {i}@{e['release']}\nPath: {r['path']}\nStatus: {r['status']}\nVerified: {'yes' if r['verified'] else 'no'}\nFiles checked: {len(r['checked_files'])}")
    return OK if r["verified"] else CONFLICT

def cmd_install(cat:dict[str,Any],root:Path,a:argparse.Namespace)->int:
    i,e=get_entry(cat,a.database_id);r=install(root,i,e)
    if a.json:emit(r)
    else:print(f"Action: {r['action']}\nStatus: {r['installation']['status']}\nPath: {r['installation']['path']}")
    return OK

def parser()->argparse.ArgumentParser:
    pr=project_root();p=argparse.ArgumentParser(description="Manage catalogued metagenomic databases.")
    p.add_argument("--catalog",default=str(pr/"config"/"database_catalog.json"));p.add_argument("--root",default=str(pr/"db"));p.add_argument("--version",action="version",version=f"%(prog)s {VERSION}")
    sp=p.add_subparsers(dest="command",required=True)
    for name,fn,helptext in (("list",cmd_list,"List catalog entries."),("info",cmd_info,"Show catalog entry."),("status",cmd_status,"Show installation status."),("verify",cmd_verify,"Verify installation."),("install",cmd_install,"Download and install database.")):
        q=sp.add_parser(name,help=helptext)
        if name!="list":q.add_argument("database_id",nargs="?")
        q.add_argument("--json",action="store_true");q.set_defaults(handler=fn)
    return p

def main(argv:list[str]|None=None)->int:
    a=parser().parse_args(argv)
    try:cat,root=context(a);return int(a.handler(cat,root,a))
    except KeyError as ex:eprint(f"ERROR: Unknown database: {ex.args[0]}");return UNKNOWN_DB
    except NoSpaceError as ex:eprint(f"ERROR: {ex}");return NO_SPACE
    except ConflictError as ex:eprint(f"ERROR: {ex}");return CONFLICT
    except ManagedError as ex:eprint(f"ERROR: {ex}");return INSTALL_ERROR if a.command=="install" else USAGE_ERROR
    except KeyboardInterrupt:eprint("ERROR: interrupted by user");return 130
    except OSError as ex:eprint(f"ERROR: {ex}");return INSTALL_ERROR
if __name__=="__main__":raise SystemExit(main())
