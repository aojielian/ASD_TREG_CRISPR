"""Optional comparison from separately obtained frozen public expression inputs.

The matrices are deliberately not distributed. This command checks exact input
identity and recomputes CPM, library consistency, and HPA values using released IDs.
"""
import argparse
import hashlib
from pathlib import Path

import numpy as np
import pandas as pd


def verify(path, expected):
    with path.open("rb") as f:
        assert hashlib.file_digest(f, "sha256").hexdigest() == expected, path.name


def run_fields(runs):
    return [f"DRR{500000+i}" for i in runs]


def expression_projection(background, count, source_index, label, group1, group2=(), minimum=3):
    counts = pd.read_csv(count, sep="\t", index_col=source_index)
    selected = [*group1, *group2]
    values = counts[selected].astype(float)
    assert values.index.is_unique and np.isfinite(values).all().all() and (values >= 0).all().all()
    totals = values.sum()
    assert (totals > 0).all()
    cpm = values.div(totals) * 1e6
    derived = pd.DataFrame({f"{label}_mean_cpm":cpm.mean(axis=1),
                            f"{label}_mean_log_expression":np.log2(cpm+1).mean(axis=1)})
    for cutoff in (.5,1,2):
        n1 = cpm[group1].ge(cutoff).sum(axis=1)
        derived[f"{label}_group1_n_ge_{cutoff:g}"] = n1
        detectable = n1.ge(minimum)
        if group2:
            n2 = cpm[list(group2)].ge(cutoff).sum(axis=1)
            derived[f"{label}_group2_n_ge_{cutoff:g}"] = n2
            detectable &= n2.ge(minimum)
        derived[f"{label}_detect_{cutoff:g}"] = detectable
    compare_projection(background, f"{label}_reference_raw_id", derived)
    print(f"{label}: verified {len(selected)} libraries and {len(values)} source genes; all released summaries agree.")


def compare_projection(background, key, derived):
    measured = background[[key, *derived.columns]].dropna(subset=[key]).drop_duplicates()
    assert measured[key].isin(derived.index).all()
    exact = derived.loc[measured[key]].reset_index(drop=True)
    measured = measured.reset_index(drop=True)
    for col in derived:
        if col.endswith(("detect_0.5","detect_1","detect_2")):
            assert measured[col].map(str).str.lower().eq(exact[col].map(str).str.lower()).all(), col
        else:
            np.testing.assert_allclose(measured[col].astype(float), exact[col].astype(float),
                                       rtol=1e-12,atol=1e-12,equal_nan=True,err_msg=col)


def check_hpa(background,path):
    hpa = pd.read_csv(path, sep="\t", usecols=["Gene","Immune cell","nTPM"])
    wanted = ["T-reg","memory CD4 T-cell","naive CD4 T-cell"]
    pivot = hpa[hpa["Immune cell"].isin(wanted)].pivot(index="Gene",columns="Immune cell",values="nTPM")
    values = pivot[wanted].add_prefix("hpa_")
    values["hpa_tconv_max_nTPM"] = pivot[wanted[1:]].max(axis=1,skipna=False)
    compare_projection(background,"hpa_reference_raw_id", values)
    print("HPA: verified T-reg and both conventional CD4 source fields and the frozen combined rule.")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("data_root",type=Path)
    p.add_argument("--itreg-counts",type=Path,required=True)
    p.add_argument("--peripheral-counts",type=Path,required=True)
    p.add_argument("--hpa",type=Path)
    args=p.parse_args()
    data=pd.read_csv(args.data_root/"data/r7/inputs/R7_GENE_BACKGROUND.tsv.gz",sep="\t",low_memory=False)
    verify(args.itreg_counts,"f240793e28c158dd820db0da123b364dc8b95841793b48ab875f0ba93f6b344a")
    verify(args.peripheral_counts,"d4e46468f9fe7c477b85106998e1d685735fcfaeb4b4db32edfbc0d08ad99f27")
    expression_projection(data,args.itreg_counts,"Gene","itreg",run_fields([405,409,413]),run_fields([406,410,414]),2)
    expression_projection(data,args.peripheral_counts,"gene_id","htreg",run_fields([418,420,422,424]))
    expression_projection(data,args.peripheral_counts,"gene_id","htconv",run_fields([417,419,421,423]))
    if args.hpa:
        verify(args.hpa,"206f4d2be7bb58731dfead00185081532b9b022164ce28e7878171e34372ccba")
        check_hpa(data,args.hpa)
    else:
        print("HPA original-file test not requested; its released model inputs are still available.")


if __name__ == "__main__":
    main()
