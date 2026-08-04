#!/usr/bin/env python3
"""
03_baseline_check.py

Control analysis. The Figure 3A metric (Spearman rho between observed post-FMT
abundance and pre_abundance * 2^predicted_log2FC) is evaluated against two
baselines that carry no model information:

  A. no-change baseline    : predict post = pre  (log2FC := 0)
  B. shuffled-FC null      : keep pre_abundance, permute the predicted log2FC

If the model does not beat these, the abundance-level rho is being driven by
pre_abundance autocorrelation rather than by learned engraftment dynamics.
Run on both the IBS (Figure 3) and Li/FAT tables so the two are comparable.

Outputs are written to Figure2/data/supp_5/.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_DIR = os.path.dirname(BASE)
OUT = os.path.join(BASE, "data", "supp_5")

FEATURES = ["average_distance", "diversity_ratio", "met_independence", "log10_pre_abund"]
SEED = 123
RF_KW = dict(n_estimators=100, max_features=3, min_samples_leaf=2,
             random_state=SEED, n_jobs=-1)


def clean(df):
    d = df.dropna(subset=FEATURES + ["log_fold_change", "pre_abundance", "post_abundance"])
    return d[np.isfinite(d.log_fold_change) & np.isfinite(d.log10_pre_abund)].copy()


def load_ibs():
    f3 = os.path.join(REPO_DIR, "Figure3")
    d = pd.read_csv(os.path.join(f3, "processed_data", "nearest_neighbors_t0.csv"))
    d[["donor", "recipient"]] = d.comparison.str.split(" vs ", expand=True)
    div = pd.read_csv(os.path.join(f3, "processed_data", "shannonDiversity.csv")).dropna()
    sc = pd.read_csv(os.path.join(f3, "processed_data", "genome_scores.tsv"), sep="\t")
    st = pd.read_csv(os.path.join(f3, "ref", "filtered_strain_taxa.tsv"), sep="\t")
    ind = (sc.merge(st, left_on="genome_name", right_on="MicrobeID")
           .rename(columns={"genome_score": "met_independence", "NCBI Taxonomy ID": "taxon"})
           .groupby("taxon", as_index=False)["met_independence"].mean())
    ab = pd.read_csv(os.path.join(f3, "processed_data", "pre_post_abundances.csv"))
    x = (d.drop(columns=["sample_id"])
         .merge(div, on=["donor", "recipient"], how="left")
         .merge(ind, on="taxon", how="left")
         .merge(ab, on=["donor", "recipient", "taxon"], how="left"))
    x["log10_pre_abund"] = np.log10(x.pre_abundance)
    return clean(x)


def loo(df):
    df = df.copy()
    df["pred_fc"] = np.nan
    for rec in df.recipient.unique():
        m = RandomForestRegressor(**RF_KW)
        m.fit(df.loc[df.recipient != rec, FEATURES], df.loc[df.recipient != rec, "log_fold_change"])
        df.loc[df.recipient == rec, "pred_fc"] = m.predict(df.loc[df.recipient == rec, FEATURES])
    df["pred_abundance"] = df.pre_abundance * 2 ** df.pred_fc
    return df


def rho_by_rec(df, pred):
    return np.array([spearmanr(g.post_abundance, g[pred]).statistic
                     for _, g in df.groupby("recipient")])


def report(name, df):
    rng = np.random.default_rng(SEED)
    model = rho_by_rec(df, "pred_abundance")
    nochange = rho_by_rec(df, "pre_abundance")
    null = []
    for _ in range(200):
        t = df.copy()
        t["shuf"] = t.pre_abundance * 2 ** rng.permutation(t.pred_fc.values)
        null.append(np.median(rho_by_rec(t, "shuf")))
    null = np.array(null)
    fc = np.array([spearmanr(g.log_fold_change, g.pred_fc).statistic
                   for _, g in df.groupby("recipient")])
    rows = [
        dict(dataset=name, metric="rho abundance, model", median=np.median(model), sd=model.std()),
        dict(dataset=name, metric="rho abundance, no-change baseline",
             median=np.median(nochange), sd=nochange.std()),
        dict(dataset=name, metric="rho abundance, shuffled-FC null",
             median=np.median(null), sd=null.std()),
        dict(dataset=name, metric="rho fold-change, model", median=np.median(fc), sd=fc.std()),
        dict(dataset=name, metric="delta vs no-change",
             median=np.median(model) - np.median(nochange), sd=np.nan),
    ]
    return pd.DataFrame(rows)


ibs = loo(load_ibs())
li = loo(clean(pd.read_csv(os.path.join(OUT, "AllMetabolicFeatures.csv"))
               .query("timepoint_day == 42")))

res = pd.concat([report("IBS (Figure 3)", ibs), report("Li/FAT day42", li)], ignore_index=True)
print(res.round(3).to_string(index=False))
res.to_csv(os.path.join(OUT, "baseline_control_IBS_vs_Li.csv"), index=False)

# how much of the abundance-level rho is just pre-abundance?
for name, d in [("IBS", ibs), ("Li", li)]:
    r = np.median([spearmanr(g.pre_abundance, g.post_abundance).statistic
                   for _, g in d.groupby("recipient")])
    v = np.median([spearmanr(g.log_fold_change, g.log10_pre_abund).statistic
                   for _, g in d.groupby("recipient")])
    print(f"{name}: median rho(pre, post) = {r:.3f}; "
          f"median rho(log2FC, log10 pre) = {v:.3f}")
