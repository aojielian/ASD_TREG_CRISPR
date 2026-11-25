"""Recompute frozen genome-wide screen margins from released inputs."""

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact
from scipy.stats.contingency import odds_ratio

PANEL = "CREBBP CUL3 DYRK1A EP300 FOXP1 KDM3B KDM6B KMT5B PHIP SIN3A SMARCA4 TBL1XR1 TLE3 ZMYND8".split()
FROZEN_FILE = {
    "chen_genomewide_analysis_ready.tsv": "Chen induced Treg FOXP3",
    "umhoefer_primary_treg_analysis_ready.tsv": "Umhoefer primary Treg FOXP3",
    "umhoefer_primary_tconv_analysis_ready.tsv": "Umhoefer primary stimulated Tconv FOXP3",
    "shifrut_asd_enrichment.tsv": "Shifrut primary CD8 TCR proliferation",
}


def write(rows, path):
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, na_rep="NA")


def binary(series):
    result = series.astype("string").str.upper().map({"1": True, "0": False, "TRUE": True, "FALSE": False})
    assert result.notna().all()
    return result.astype(bool)


def baseline(data, shifrut, out):
    processed = data / "data/processed"
    rows = []
    margins = [(14, 228, 243, 18481), (5, 28, 26, 1290), (4, 29, 34, 1282), (11, 231, 965, 17887)]
    for (name, label), reference in zip(FROZEN_FILE.items(), margins):
        frame = pd.read_csv(processed / name, sep="\t", low_memory=False)
        if "sfari_score1" in frame:
            case, hit = binary(frame.sfari_score1), binary(frame.screen_hit)
            assert not frame.gene.duplicated().any()
            counts = ((case & hit).sum(), (case & ~hit).sum(), (~case & hit).sum(), (~case & ~hit).sum())
            if name.startswith("chen"):
                assert sorted(frame.loc[case & hit, "gene"].tolist()) == PANEL
        else:
            old = shifrut
            genes = pd.read_csv(old, sep="\t")
            assert not genes.approved_symbol.duplicated().any()
            case, hit = binary(genes.sfari_score1), binary(genes.SHIFRUT_HIT)
            author_rule = genes.neg_fdr.le(.05) | genes.pos_fdr.le(.05)
            assert hit.eq(author_rule).all()
            counts = ((case & hit).sum(), (case & ~hit).sum(), (~case & hit).sum(), (~case & ~hit).sum())
            source = frame.iloc[0]
            assert counts == tuple(int(source[k]) for k in ["asd_hits", "asd_nonhits", "other_hits", "other_nonhits"])
            name = "data/r10/reference/" + old.name + "; " + name
        if counts != reference:
            raise ValueError(f"{label} {counts}, frozen {reference}")
        table = np.array(counts).reshape(2, 2)
        est = odds_ratio(table, kind="conditional")
        low, high = est.confidence_interval()
        one = fisher_exact(table, alternative="greater").pvalue
        two = fisher_exact(table, alternative="two-sided").pvalue
        rows.append(dict(screen=label, N=sum(counts), asd_hit=counts[0], asd_nonhit=counts[1],
                         other_hit=counts[2], other_nonhit=counts[3], conditional_OR=est.statistic,
                         exact_CI95_low=low, exact_CI95_high=high, p_one_sided_greater=one,
                         p_two_sided=two, source=name))
    write(rows, out / "R10_R9_BASELINE_RECHECK.tsv")
