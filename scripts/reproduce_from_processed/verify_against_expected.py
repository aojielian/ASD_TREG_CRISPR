#!/usr/bin/env python3
from __future__ import annotations
import csv, math, sys
def read(path):
    with open(path,encoding="utf-8",newline="") as h: return list(csv.DictReader(h,delimiter="\t"))
def main():
    expected=read(sys.argv[1]); observed={}
    for path in sys.argv[2:-1]:
        for r in read(path): observed[r["result_id"]]=r["observed_value"]
    rows=[]
    for r in expected:
        key=r["result_id"]; exp=r["expected_value"]; obs=observed.get(key); tol=float(r["tolerance"] or 0); ok=False
        if obs is not None:
            if r["value_type"] in {"integer","string"}: ok=obs==exp
            else:
                try: ok=abs(float(obs)-float(exp))<=tol
                except ValueError: ok=False
        delta=""
        if obs is not None and r["value_type"]=="numeric":
            try: delta=str(float(obs)-float(exp))
            except ValueError: delta="NA"
        rows.append({"result_id":key,"expected_value":exp,"observed_value":obs or "","tolerance":r["tolerance"],"delta":delta,"status":"PASS" if ok else "FAIL","source_basis":r["source_basis"]})
    with open(sys.argv[-1],"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=list(rows[0]),delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
    status="PASS" if all(r["status"]=="PASS" for r in rows) else "FAIL"
    print("EXPECTED_RESULTS_VERIFICATION="+status)
    if status!="PASS": sys.exit(1)
if __name__=="__main__": main()
