#!/usr/bin/env python3
"""
04_diagnostics.py

Why does the model transfer at the abundance level but not at the fold-change
level in Li/FAT? Three diagnostics:

  1. Feature ablation on both datasets. Single-feature (log10_pre_abund) model
     vs full model, LOO fold-change rho. Separates "regression to the mean on
     pre-abundance" from "metabolic features carry engraftment information".
  2. Allogenic-only vs autologous-only Li models across all four timepoints.
  3. Dataset comparability: feature distributions and perturbation magnitude.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor
import importlib.util

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "data", "supp_5")
spec = importlib.util.spec_from_file_location(
    "bc", os.path.join(BASE, "scripts", "03_baseline_check.py"))

FEATURES = ["average_distance", "diversity_ratio", "met_independence", "log10_pre_abund"]
RF_KW = dict(n_estimators=100, min_samples_leaf=2, random_state=123, n_jobs=-1)


def clean(df):
    d = df.dropna(subset=FEATURES + ["log_fold_change", "pre_abundance", "post_abundance"])
    return d[np.isfinite(d.log_fold_change) & np.isfinite(d.log10_pre_abund)].copy()


def loo_fc_rho(df, feats):
    df = df.copy()
    df["p"] = np.nan
    for rec in df.recipient.unique():
        tr, te = df[df.recipient != rec], df[df.recipient == rec]
        m = RandomForestRegressor(**{**RF_KW, "max_features": min(3, len(feats))})
        m.fit(tr[feats], tr.log_fold_change)
        df.loc[df.recipient == rec, "p"] = m.predict(te[feats])
    per = np.array([spearmanr(g.log_fold_change, g.p).statistic for _, g in df.groupby("recipient")])
    abd = np.array([spearmanr(g.post_abundance, g.pre_abundance * 2 ** g.p).statistic
                    for _, g in df.groupby("recipient")])
    return np.median(per), per.std(), np.median(abd)


# reuse the IBS loader from 03
mod = {"__file__": os.path.join(BASE, "scripts", "03_baseline_check.py")}
exec(open(os.path.join(BASE, "scripts", "03_baseline_check.py")).read().split("ibs = loo(")[0],
     mod)
ibs = mod["load_ibs"]()
li_all = pd.read_csv(os.path.join(OUT, "AllMetabolicFeatures.csv"))
li = clean(li_all.query("timepoint_day == 42"))

print("=== 1. feature ablation, LOO fold-change rho ===")
rows = []
for name, d in [("IBS", ibs), ("Li day42", li)]:
    for label, feats in [("full (4 features)", FEATURES),
                         ("log10_pre_abund only", ["log10_pre_abund"]),
                         ("metabolic only (no pre-abund)",
                          ["average_distance", "diversity_ratio", "met_independence"])]:
        med, sd, abd = loo_fc_rho(d, feats)
        rows.append(dict(dataset=name, features=label, rho_fc=med, sd_fc=sd, rho_abundance=abd))
abl = pd.DataFrame(rows)
print(abl.round(3).to_string(index=False))
abl.to_csv(os.path.join(OUT, "feature_ablation.csv"), index=False)

print("\n=== 2. Li by fmt_type and timepoint, LOO fold-change rho ===")
rows = []
for day in sorted(li_all.timepoint_day.unique()):
    d = clean(li_all[li_all.timepoint_day == day])
    for grp in ["Allogenic", "Autologous", "all"]:
        sub = d if grp == "all" else d[d.fmt_type == grp]
        med, sd, abd = loo_fc_rho(sub, FEATURES)
        rows.append(dict(day=day, group=grp, n_recipients=sub.recipient.nunique(),
                         rho_fc=med, sd_fc=sd, rho_abundance=abd))
grp_df = pd.DataFrame(rows)
print(grp_df.round(3).to_string(index=False))
grp_df.to_csv(os.path.join(OUT, "Li_by_group_and_timepoint.csv"), index=False)

print("\n=== 3. dataset comparability ===")
comp = []
for name, d in [("IBS", ibs), ("Li day42", li)]:
    comp.append(dict(dataset=name,
                     n_rows=len(d), n_recipients=d.recipient.nunique(), n_taxa=d.taxon.nunique(),
                     sd_log2FC=d.log_fold_change.std(),
                     iqr_log2FC=d.log_fold_change.quantile(.75) - d.log_fold_change.quantile(.25),
                     med_abs_log2FC=d.log_fold_change.abs().median(),
                     mean_avg_distance=d.average_distance.mean(),
                     mean_div_ratio=d.diversity_ratio.mean(),
                     mean_met_indep=d.met_independence.mean(),
                     mean_log10_pre=d.log10_pre_abund.mean()))
comp = pd.DataFrame(comp)
print(comp.round(3).to_string(index=False))
comp.to_csv(os.path.join(OUT, "dataset_comparability.csv"), index=False)

# how much perturbation is there at all?
print("\nfraction of taxa with |log2FC| > 1:")
for name, d in [("IBS", ibs), ("Li day42", li)]:
    print(f"  {name}: {(d.log_fold_change.abs() > 1).mean():.3f}")
for day in sorted(li_all.timepoint_day.unique()):
    d = clean(li_all[li_all.timepoint_day == day])
    for g in ["Allogenic", "Autologous"]:
        s = d[d.fmt_type == g]
        print(f"  Li day{day} {g}: {(s.log_fold_change.abs() > 1).mean():.3f} "
              f"(sd log2FC {s.log_fold_change.std():.2f})")
