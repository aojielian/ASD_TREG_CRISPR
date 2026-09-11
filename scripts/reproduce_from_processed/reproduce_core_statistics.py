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
def conditional_or(a,b_,c,d):
    N=a+b_+c+d; m=a+c; nn=b_+d; draws=a+b_; lo=max(0,draws-nn); hi=min(draws,m)
    if a==hi: return float("inf")
    if a==lo: return 0.0
    support=list(range(lo,hi+1)); logdc=[]
    for k in support:
        logdc.append(math.lgamma(m+1)-math.lgamma(k+1)-math.lgamma(m-k+1)+math.lgamma(nn+1)-math.lgamma(draws-k+1)-math.lgamma(nn-draws+k+1))
    def mean(ncp):
        z=[x+math.log(ncp)*k for x,k in zip(logdc,support)]; mx=max(z); w=[math.exp(x-mx) for x in z]
        return sum(k*x for k,x in zip(support,w))/sum(w)
    def brent(f,left,right,tol=2**-13):
        fa=f(left); fb=f(right); c=right; fc=fb; step=right-left; old=step
        for _ in range(1000):
            if (fb>0 and fc>0) or (fb<0 and fc<0): c=left; fc=fa; step=old=right-left
            if abs(fc)<abs(fb): left,right,c=right,c,right; fa,fb,fc=fb,fc,fb
            tol1=2*2.220446049250313e-16*abs(right)+tol/2; xm=(c-right)/2
            if abs(xm)<=tol1 or fb==0: return right
            if abs(old)>=tol1 and abs(fa)>abs(fb):
                s=fb/fa
                if left==c: pp=2*xm*s; qq=1-s
                else:
                    qq=fa/fc; rr=fb/fc; pp=s*(2*xm*qq*(qq-rr)-(right-left)*(rr-1)); qq=(qq-1)*(rr-1)*(s-1)
                if pp>0: qq=-qq
                else: pp=-pp
                if 2*pp<min(3*xm*qq-abs(tol1*qq),abs(old*qq)): old=step; step=pp/qq
                else: step=xm; old=step
            else: step=xm; old=step
            left=right; fa=fb; right += step if abs(step)>tol1 else math.copysign(tol1,xm); fb=f(right)
        return right
    return 1/brent(lambda t:mean(1/t)-a,2.220446049250313e-16,1)
def fisher(a,b_,c,d):
    N=a+b_+c+d; r=a+b_; col=a+c; den=math.comb(N,col); p=0.0
    for k in range(a,min(r,col)+1): p += math.comb(r,k)*math.comb(N-r,col-k)/den
    return conditional_or(a,b_,c,d),p
def fisher_two(a,b_,c,d):
    N=a+b_+c+d; r=a+b_; col=a+c; den=math.comb(N,col)
    probs=[(k,math.comb(r,k)*math.comb(N-r,col-k)/den) for k in range(max(0,col-(N-r)),min(r,col)+1)]
    obs=next(p for k,p in probs if k==a)
    return sum(p for k,p in probs if p<=obs+1e-15)
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
    return .5*math.erfc((W-mu-cc)/math.sqrt(2*var))
def hyper(N,H,A):
    den=math.comb(N,A); lo=max(0,A-(N-H)); hi=min(A,H)
    return {k:math.comb(H,k)*math.comb(N-H,A-k)/den for k in range(lo,hi+1)}
def convolution(rows,idcol):
    groups=defaultdict(list)
    for r in rows: groups[r[idcol]].append(r)
    pm={0:1.0}; obs=0; exp=0.0
    for rs in groups.values():
        N=len(rs); H=sum(b(r["screen_hit"]) for r in rs); A=sum(b(r["sfari_score1"]) for r in rs); K=sum(b(r["sfari_score1"]) and b(r["screen_hit"]) for r in rs)
        obs+=K; exp+=A*H/N; q=hyper(N,H,A); nxt=defaultdict(float)
        for x,px in pm.items():
            for k,pk in q.items(): nxt[x+k]+=px*pk
        pm=nxt
    return obs,exp,sum(p for k,p in pm.items() if k>=obs)
def holm_two(pvalues):
    """Return direct Holm-adjusted values for exactly two named tests."""
    ordered = sorted(pvalues.items(), key=lambda item: (item[1], item[0]))
    first = min(1.0, 2.0 * ordered[0][1])
    second = max(first, ordered[1][1])
    return {ordered[0][0]: first, ordered[1][0]: second}
def write(rows,path):
    with open(path,"w",encoding="utf-8",newline="") as h:
        w=csv.DictWriter(h,fieldnames=["result_id","observed_value"],delimiter="\t",lineterminator="\n"); w.writeheader(); w.writerows(rows)
def main():
    snap=Path(sys.argv[1]); p=snap/"data/processed"; out=[]

    chen=read(p/"chen_genomewide_analysis_ready.tsv")
    ca=[r for r in chen if b(r["sfari_score1"])]; co=[r for r in chen if not b(r["sfari_score1"])]
    add(out,"CHEN_FULL_UNIVERSE_ROWS",len(chen))
    a=sum(b(r["screen_hit"]) for r in ca); bb=len(ca)-a; c=sum(b(r["screen_hit"]) for r in co); d=len(co)-c
    add(out,"CHEN_PRIMARY_SFARI_HITS",a); add(out,"CHEN_PRIMARY_SFARI_NOT_HITS",bb); add(out,"CHEN_PRIMARY_OTHER_HITS",c); add(out,"CHEN_PRIMARY_OTHER_NOT_HITS",d)
    o,pv=fisher(a,bb,c,d); add(out,"CHEN_PRIMARY_OR",o); add(out,"CHEN_PRIMARY_FISHER_P",pv)
    add(out,"CHEN_SCREEN_STRENGTH_P",rank_p([n(r["screen_strength"]) for r in ca],[n(r["screen_strength"]) for r in co]))
    add(out,"CHEN_ABS_LFC_P",rank_p([n(r["abs_author_lfc"]) for r in ca],[n(r["abs_author_lfc"]) for r in co]))

    treg=read(p/"umhoefer_primary_treg_analysis_ready.tsv"); ta=[r for r in treg if b(r["sfari_score1"])]; to=[r for r in treg if not b(r["sfari_score1"])]
    a=sum(b(r["screen_hit"]) for r in ta); bb=len(ta)-a; c=sum(b(r["screen_hit"]) for r in to); d=len(to)-c
    o,pv=fisher(a,bb,c,d); add(out,"UMHOEFER_TREG_FULL_UNIVERSE_ROWS",len(treg)); add(out,"UMHOEFER_TREG_SFARI_HITS",a); add(out,"UMHOEFER_TREG_SFARI_NOT_HITS",bb); add(out,"UMHOEFER_TREG_OTHER_HITS",c); add(out,"UMHOEFER_TREG_OTHER_NOT_HITS",d); add(out,"UMHOEFER_TREG_OR",o); add(out,"UMHOEFER_TREG_FISHER_P",pv); add(out,"UMHOEFER_TREG_SCREEN_STRENGTH_P",rank_p([n(r["screen_strength"]) for r in ta],[n(r["screen_strength"]) for r in to]))
    tc=read(p/"umhoefer_primary_tconv_analysis_ready.tsv"); tca=[r for r in tc if b(r["sfari_score1"])]; tco=[r for r in tc if not b(r["sfari_score1"])]
    a=sum(b(r["screen_hit"]) for r in tca); bb=len(tca)-a; c=sum(b(r["screen_hit"]) for r in tco); d=len(tco)-c; o,pv=fisher(a,bb,c,d); add(out,"UMHOEFER_TCONV_SFARI_HITS",a); add(out,"UMHOEFER_TCONV_OTHER_HITS",c); add(out,"UMHOEFER_TCONV_OR",o); add(out,"UMHOEFER_TCONV_FISHER_P",pv)

    cs=read(p/"cross_screen_concordance.tsv"); genes=[r for r in cs if r.get("row_type")=="COMMON_GENE"]; x=sum(b(r["chen_hit"]) and b(r["umhoefer_treg_hit"]) for r in genes); y=sum(b(r["chen_hit"]) and not b(r["umhoefer_treg_hit"]) for r in genes); z=sum(not b(r["chen_hit"]) and b(r["umhoefer_treg_hit"]) for r in genes); w=len(genes)-x-y-z
    add(out,"CROSS_SCREEN_COMMON_ROWS",len(genes)); add(out,"CROSS_SCREEN_BOTH_HITS",x); o,_=fisher(x,y,z,w); add(out,"CROSS_SCREEN_OR",o); add(out,"CROSS_SCREEN_P",fisher_two(x,y,z,w))
    score1=read(p/"chen_score1_screen_results.tsv"); treg_by_gene={r["gene"]:r for r in treg}
    common_score1=[r for r in score1 if b(r["chen_measured"]) and b(r["umhoefer_treg_measured"]) and r["gene"] in treg_by_gene]
    both_names=sorted(r["gene"] for r in common_score1 if r["chen_hit"] in {"YES","TRUE"} and b(treg_by_gene[r["gene"]]["screen_hit"]))
    treg_only_names=sorted(r["gene"] for r in common_score1 if r["chen_hit"] not in {"YES","TRUE"} and b(treg_by_gene[r["gene"]]["screen_hit"]))
    neither_names=sorted(r["gene"] for r in common_score1 if r["chen_hit"] not in {"YES","TRUE"} and not b(treg_by_gene[r["gene"]]["screen_hit"]))
    add(out,"CROSS_SCREEN_ASD_SCORE1_COMMON_ROWS",len(common_score1)); add(out,"CROSS_SCREEN_ASD_BOTH_HIT_COUNT",len(both_names)); add(out,"CROSS_SCREEN_ASD_BOTH_HIT_NAMES",";".join(both_names)); add(out,"CROSS_SCREEN_ASD_TREG_ONLY_COUNT",len(treg_only_names)); add(out,"CROSS_SCREEN_ASD_TREG_ONLY_NAMES",";".join(treg_only_names)); add(out,"CROSS_SCREEN_ASD_NEITHER_COUNT",len(neither_names))
    rem=read(p/"observed_hit_removal_sensitivity.tsv")
    for r,prefix in [(next(r for r in rem if r["screen"]=="CHEN"),"FOXP1_ZMYND8_CHEN"),(next(r for r in rem if r["screen"]=="UMHOEFER_TREG"),"FOXP1_ZMYND8_TREG")]:
        add(out,prefix+"_REMAINING_HITS",r["K"]); add(out,prefix+"_P",r["one_sided_fisher_p"])

    sh=read(p/"chen_shifrut_common_universe.tsv"); add(out,"CHEN_SHIFRUT_COMMON_ROWS",len(sh))
    for cat in ["FOXP3_ONLY","CONCORDANT","TCR_ONLY"]: add(out,"CHEN_SHIFRUT_"+cat,sum(r["derived_category"]==cat for r in sh))
    for cat in ["FOXP3_ONLY","CONCORDANT","TCR_ONLY"]: add(out,"CHEN_SHIFRUT_ASD_"+cat,sum(r["derived_category"]==cat and b(r["SFARI_SCORE1"]) for r in sh))
    A=sum(b(r["SFARI_SCORE1"]) for r in sh); S=sum(int(r["D"]) for r in sh if b(r["SFARI_SCORE1"])); add(out,"DIRECT_PREFERENCE_S_OBSERVED",S); N=len(sh); np=sum(int(r["D"])==1 for r in sh); nz=sum(int(r["D"])==0 for r in sh); nm=sum(int(r["D"])==-1 for r in sh); den=math.comb(N,A); tail=0
    for kp in range(min(np,A)+1):
        for km in range(min(nm,A-kp)+1):
            kz=A-kp-km
            if 0<=kz<=nz:
                ways=math.comb(np,kp)*math.comb(nm,km)*math.comb(nz,kz)
                if kp-km>=S: tail+=ways
    add(out,"DIRECT_PREFERENCE_EXACT_P",tail/den)

    # Class conditioning uses the complete full universes; joint conditioning uses complete cases.
    cclass=[]
    for r in chen:
        z=dict(r); z["final_class_stratum_id"]=r["class_label"]+"|"+("INSIDE" if r["library_member_umhoefer_treg"]=="TRUE" else "OUTSIDE"); cclass.append(z)
    tclass=[]
    for r in treg:
        z=dict(r); z["final_class_stratum_id"]=r["class_label"]; tclass.append(z)
    class_values = {}
    for rows,prefix,screen in [(cclass,"CHEN_CLASS","CHEN"),(tclass,"TREG_CLASS","UMHOEFER_TREG")]:
        o,e,pv=convolution(rows,"final_class_stratum_id"); class_values[screen]=(prefix,o,e,pv)
    class_holm = holm_two({screen: values[3] for screen, values in class_values.items()})
    for screen,(prefix,o,e,pv) in class_values.items():
        add(out,prefix+"_OBSERVED_K",o); add(out,prefix+"_EXPECTED_K",e); add(out,prefix+"_OE",o/e); add(out,prefix+"_EXACT_P",pv); add(out,prefix+"_HOLM_P",class_holm[screen])
    if sorted(both_names) != ["FOXP1", "ZMYND8"] or sorted(treg_only_names) != ["CTCF", "PTEN", "SATB1"]:
        raise ValueError("common ASD score-1 gene partition does not match the released result")
    for fn,prefix,screen in [("chen_joint_conditioning_complete_case.tsv","CHEN_JOINT","CHEN"),("umhoefer_treg_joint_conditioning_complete_case.tsv","TREG_JOINT","UMHOEFER_TREG")]:
        rows=read(p/fn); o,e,pv=convolution(rows,"final_joint_stratum_id"); class_values[screen]=(prefix,o,e,pv)
    joint_holm = holm_two({screen: class_values[screen][3] for screen in ("CHEN", "UMHOEFER_TREG")})
    for screen in ("CHEN", "UMHOEFER_TREG"):
        prefix,o,e,pv=class_values[screen]; add(out,prefix+"_OBSERVED_K",o); add(out,prefix+"_EXPECTED_K",e); add(out,prefix+"_OE",o/e); add(out,prefix+"_EXACT_P",pv); add(out,prefix+"_HOLM_P",joint_holm[screen])
    for fn,prefix,full_fn in [("chen_joint_conditioning_complete_case.tsv","CHEN_JOINT","chen_genomewide_analysis_ready.tsv"),("umhoefer_treg_joint_conditioning_complete_case.tsv","TREG_JOINT","umhoefer_primary_treg_analysis_ready.tsv")]:
        rows=read(p/fn); full=read(p/full_fn); complete_asd=sum(b(r["sfari_score1"]) for r in rows); original_asd=sum(b(r["sfari_score1"]) for r in full)
        add(out,prefix+"_COMPLETE_CASE_ASD_MEASURED",complete_asd); add(out,prefix+"_ORIGINAL_ASD_MEASURED",original_asd); add(out,prefix+"_ASD_COVERAGE",complete_asd/original_asd)
    add(out,"CHEN_CLASS_COMPLETE_CASE_ROWS",len(read(p/"chen_joint_conditioning_complete_case.tsv"))); add(out,"TREG_CLASS_COMPLETE_CASE_ROWS",len(read(p/"umhoefer_treg_joint_conditioning_complete_case.tsv")))
    for fn,prefix in [("chen_joint_conditioning_complete_case.tsv","CHEN_COMPLETE_CASE"),("umhoefer_treg_joint_conditioning_complete_case.tsv","TREG_COMPLETE_CASE")]:
        rows=read(p/fn); A=sum(b(r["sfari_score1"]) for r in rows); O=len(rows)-A; ah=sum(b(r["sfari_score1"]) and b(r["screen_hit"]) for r in rows); oh=sum(not b(r["sfari_score1"]) and b(r["screen_hit"]) for r in rows); o,pv=fisher(ah,A-ah,oh,O-oh); add(out,prefix+"_OR",o); add(out,prefix+"_P",pv)

    add(out,"M1_M10_ACCEPTED_ROWS",len(read(p/"chen_M1_M10_canonical_memberships.tsv"))); add(out,"M1_M10_UNRESOLVED_ROWS",len(read(p/"M1_M10_unresolved_source_entries.tsv")))
    loo=read(p/"loo_mouse_analysis_ready.tsv"); A=sum(b(r["sfari_score1"]) for r in loo); O=len(loo)-A; ah=sum(b(r["sfari_score1"]) and b(r["LOO_FOXP3_HIT"]) for r in loo); oh=sum(not b(r["sfari_score1"]) and b(r["LOO_FOXP3_HIT"]) for r in loo); o,pv=fisher(ah,A-ah,oh,O-oh); released=read(p/"loo_primary_enrichment.tsv")[0]; add(out,"LOO_ANALYSIS_READY_ROWS",len(loo)); add(out,"LOO_ASD_HITS",ah); add(out,"LOO_OTHER_HITS",oh); add(out,"LOO_REFERENCE_R_CONDITIONAL_OR",released["odds_ratio"]); add(out,"LOO_FISHER_P",pv)
    write(out,sys.argv[2])
if __name__=="__main__": main()
