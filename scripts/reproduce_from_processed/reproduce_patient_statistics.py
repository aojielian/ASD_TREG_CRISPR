#!/usr/bin/env python3
from __future__ import annotations
import csv, math, sys
from collections import defaultdict
from pathlib import Path

def read(path):
    with open(path,encoding="utf-8",newline="") as h: return list(csv.DictReader(h,delimiter="\t"))
def n(v):
    try: return float(v) if v not in (None,"","NA") else None
    except ValueError: return None
def b(v): return str(v).upper() in {"TRUE","YES","1"}
def add(out,key,value): out.append({"result_id":key,"observed_value":str(value)})
def rank_p(x,y):
    z=sorted([(v,0) for v in x]+[(v,1) for v in y]); ranks=[0.0]*len(z); ties=defaultdict(int); i=0
    while i<len(z):
        j=i
        while j+1<len(z) and z[j+1][0]==z[i][0]: j+=1
        rr=(i+1+j+1)/2
        for k in range(i,j+1): ranks[k]=rr
        ties[z[i][0]]=j-i+1; i=j+1
    n1=len(x); n2=len(y); N=n1+n2; W=sum(r for r,(_,g) in zip(ranks,z) if g==0); mu=n1*(N+1)/2
    var=n1*n2/12*((N+1)-sum(t**3-t for t in ties.values())/(N*(N-1)))
    cc=.5 if W>mu else -.5 if W<mu else 0
    return math.erfc(abs(W-mu-cc)/math.sqrt(2*var))
def signflip(d):
    obs=sum(d); ge=0
    for mask in range(1<<len(d)):
        s=sum((1 if mask>>i&1 else -1)*x for i,x in enumerate(d))
        if abs(s)>=abs(obs)-1e-12: ge+=1
    return ge/(1<<len(d))
def write(rows,path):
    with open(path,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=["result_id","observed_value"],delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
def main():
    snap=Path(sys.argv[1]); p=snap/"data/processed"; out=[]
    wb=read(p/"whole_blood_pair_program_differences.tsv"); add(out,"WHOLE_BLOOD_PAIR_PROGRAM_ROWS",len(wb))
    programs=sorted({r["program_id"] for r in wb},key=lambda x:int(x.split()[-1]))
    raw=[signflip([n(r["ASD_minus_sibling"]) for r in wb if r["program_id"]==q]) for q in programs]
    order=sorted(range(len(raw)),key=lambda i:raw[i]); q=[1.0]*len(raw)
    for rank,i in enumerate(order,1): q[i]=min(1.0,raw[i]*len(raw)/rank)
    for j in range(len(order)-2,-1,-1): q[order[j]]=min(q[order[j]],q[order[j+1]])
    add(out,"WHOLE_BLOOD_MIN_Q_BH10",min(q)); add(out,"WHOLE_BLOOD_Q_LT_005_COUNT",sum(x<.05 for x in q))

    mor=read(p/"moreno_gene_level_analysis_ready.tsv"); add(out,"MORENO_VALID_EFFECT_ROWS",len(mor)); add(out,"MORENO_DEG_ROWS",sum(b(r["deg_fdr05"]) for r in mor))
    ps=[]; eligible=[]; program_sizes={r["program_id"]:float(r["source_size"]) for r in read(p/"moreno_program_results.tsv")}
    for i in range(1,11):
        x=[n(r["signed_effect"]) for r in mor if b(r[f"module_{i}_member"])]
        y=[n(r["signed_effect"]) for r in mor if not b(r[f"module_{i}_member"])]
        ok=len(x)>=20 and len(x)/program_sizes[f"Module {i}"]>=.70
        ps.append(rank_p(x,y) if ok else 1.0); eligible.append(ok)
    order=sorted(range(10),key=lambda i:ps[i]); qs=[1.0]*10
    for rank,i in enumerate(order,1): qs[i]=min(1.0,ps[i]*10/rank)
    for j in range(8,-1,-1): qs[order[j]]=min(qs[order[j]],qs[order[j+1]])
    supported=[f"Module {i}" for i in range(1,11) if qs[i-1]<.05]
    add(out,"MORENO_SUPPORTED_MODULES",";".join(supported))
    for i in [1,6,8,9]: add(out,f"MORENO_M{i}_P",ps[i-1]); add(out,f"MORENO_M{i}_Q",qs[i-1])
    add(out,"MORENO_ASD14_ROWS",sum(b(r["ASD14_MEMBER"]) for r in mor)); add(out,"MORENO_ASD14_DEG",sum(b(r["ASD14_MEMBER"]) and b(r["deg_fdr05"]) for r in mor)); add(out,"MORENO_OTHER_DEG",sum((not b(r["ASD14_MEMBER"])) and b(r["deg_fdr05"]) for r in mor)); add(out,"MORENO_ASD14_FISHER_P",1)
    write(out,sys.argv[2])
if __name__=="__main__": main()
