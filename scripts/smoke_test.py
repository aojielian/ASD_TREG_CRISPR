"""Verify an unpacked v1.1.0 data root or ZIP against this checkout."""
import argparse
import csv
import hashlib
import os
from pathlib import Path, PurePosixPath
import platform
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile

REPO = Path(__file__).resolve().parents[1]
VERSION = "ASD_TREG_CRISPR_processed_data_v1.1.0"


def python_command():
    # Keep the active Intel interpreter architecture on universal macOS Python.
    cli = [sys.executable]
    if platform.system() == "Darwin" and platform.machine() == "x86_64":
        cli = ["/usr/bin/arch","-x86_64",*cli]
    return cli


def verify_archive_sidecar(path):
    sidecar=path.with_suffix(path.suffix+".sha256")
    if sidecar.exists():
        expected,name=sidecar.read_text().split(maxsplit=1)
        assert name.strip().lstrip("*")==path.name
        with path.open("rb") as f:
            assert hashlib.file_digest(f,"sha256").hexdigest()==expected


def extract(path, temporary):
    verify_archive_sidecar(path)
    with zipfile.ZipFile(path) as archive:
        assert archive.testzip() is None
        seen=set()
        for item in archive.infolist():
            rel=PurePosixPath(item.filename)
            if (rel.is_absolute() or ".." in rel.parts or rel.parts[0]!=VERSION
                or "\\" in item.filename or stat.S_ISLNK(item.external_attr >> 16)):
                raise ValueError(f"Unsafe archive layout at {item.filename}")
            if rel in seen:
                raise ValueError(f"Duplicate archive member: {item.filename}")
            seen.add(rel)
            dst=temporary/rel
            if item.is_dir():
                dst.mkdir(parents=True,exist_ok=True)
            else:
                dst.parent.mkdir(parents=True,exist_ok=True)
                with archive.open(item) as source, dst.open("xb") as target:
                    shutil.copyfileobj(source,target)
    return temporary/VERSION


def tsv(path):
    with path.open() as f:
        return list(csv.DictReader(f,delimiter="\t"))


def headline(out):
    for assertion in tsv(REPO/"expected_results.tsv"):
        table = out/"gene_property_adjustment/results"/assertion["table"]
        rows = tsv(table)
        match = dict(x.split("=",1) for x in assertion["selector"].split("|"))
        answer = [r for r in rows if all(r.get(field)==value for field,value in match.items())]
        assert len(answer)==1, assertion
        observed=answer[0][assertion["field"]]
        expected=assertion["expected"]
        if assertion["tolerance"]:
            assert abs(float(observed)-float(expected))<=float(assertion["tolerance"]), (assertion,observed)
        else:
            assert observed==expected, (assertion,observed)


def go(data,out,verify_r,rscript,figures):
    assert (data/"VERSION.txt").read_text().strip()==VERSION
    px=python_command()
    env=os.environ.copy()
    env.pop("PYTHONPATH",None)
    env.pop("PYTHONSTARTUP",None)
    env.update(PYTHONDONTWRITEBYTECODE="1",OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1",MKL_NUM_THREADS="1")
    tasks = [
        ("Release integrity", [*px,"-I",str(REPO/"tests/verify_release_manifest.py"),str(data),
                               "--report",str(out/"TSV_SCHEMA_QA.tsv")]),
        ("Legacy statistics", [*px,str(REPO/"scripts/legacy/run_legacy_reproduction.py"),str(data),
                               str(out/"legacy"),*(["--rscript",rscript] if verify_r else [])]),
        ("Gene-property model grid", [*px,str(REPO/"scripts/gene_property_adjustment/reproduce.py"),"--data-dir",str(data),
                           "--out-dir",str(out/"gene_property_adjustment"),*(["--verify-r"] if verify_r else [])]),
        ("Gene-property figure and mechanism sources", [*px,"-I",str(REPO/"scripts/experimental_evidence/verify_release_evidence.py"),
                 str(data/"data/r7"),"--report",str(out/"GENE_PROPERTY_EVIDENCE_QA.tsv")]),
    ]
    tasks.append(("Zhu compact evidence", [*px,"-I",str(REPO/"scripts/zhu/extract_zhu_evidence.py"),
                 "verify","--derived-dir",str(data/"data/r8/zhu_derived"),"--candidate-master",str(data/"data/r7/mechanism/CANDIDATE_DIRECTION_MECHANISM.tsv")]))
    if figures:
        env.update(MPLBACKEND="Agg", MPLCONFIGDIR=str(out/"matplotlib_cache"))
        tasks.extend([
            ("Main figures", [*px,str(REPO/"figures/render_main.py"),str(data),str(out/"figures")]),
            ("Supplementary figures", [*px,str(REPO/"figures/render_supp.py"),str(data),str(out/"figures")]),
        ])
    for label, cmd in tasks:
        print(f"Running {label}",flush=True)
        subprocess.run(cmd,env=env,check=True)
    headline(out)
    print("All release-only core tests and headline expectations passed.")
    if not verify_r:
        print("The optional independent base-R profiles and exact historical Loo OR check were not requested.")
    print("Original provider files and manuscript page underlays remain outside this smoke test.")


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("data",type=Path)
    p.add_argument("--out-dir",type=Path,required=True)
    p.add_argument("--verify-r",action="store_true")
    p.add_argument("--figures",action="store_true",help="render the four main panels and all six supplementary figures")
    p.add_argument("--rscript",default="Rscript",help="R executable for the optional legacy independent check")
    a=p.parse_args()
    out=a.out_dir.resolve()
    data=a.data.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("Specify an empty output directory.")
    out.mkdir(parents=True,exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="asd-treg-release-") as temp:
        root=extract(data,Path(temp)) if data.is_file() else data
        go(root,out,a.verify_r,a.rscript,a.figures)


if __name__=="__main__":
    main()
