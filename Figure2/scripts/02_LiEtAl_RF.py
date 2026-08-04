#!/usr/bin/env python3
"""
02_LiEtAl_RF.py

Two analyses, mirroring the IBS analysis used for Figure 3:

  (1) Li-native model. Random forest trained on the Li/FAT feature table with
      leave-one-recipient-out CV, same four features as Figure 3AB.
  (2) Cross-dataset transfer. The IBS model (retrained on the full IBS table,
      identical hyperparameters) applied to Li without ever seeing Li data.

Reported metric matches Figure 3A: per-recipient Spearman rho between observed
post-FMT abundance and predicted post-FMT abundance (pre * 2^predicted log2FC),
plus the more conservative fold-change-level rho.

Outputs are written to Figure2/data/supp_5/.
"""

import os
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, mannwhitneyu
from sklearn.ensemble import RandomForestRegressor

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
REPO_DIR = os.path.dirname(BASE)
OUT = os.path.join(BASE, "data", "supp_5")
os.makedirs(OUT, exist_ok=True)

FEATURES = ["average_distance", "diversity_ratio", "met_independence", "log10_pre_abund"]
PRIMARY_DAY = 84   # ~2.8 months; closest Li analogue to the IBS 6-month post-FMT sample
SEED = 123

RF_KW = dict(n_estimators=100, max_features=3, min_samples_leaf=2,
             random_state=SEED, n_jobs=-1)


def fit_rf(df, seed=SEED):
    m = RandomForestRegressor(**{**RF_KW, "random_state": seed})
    m.fit(df[FEATURES], df["log_fold_change"])
    return m


def clean(df):
    d = df.dropna(subset=FEATURES + ["log_fold_change", "pre_abundance", "post_abundance"])
    return d[np.isfinite(d.log_fold_change) & np.isfinite(d.log10_pre_abund)].copy()


def per_recipient_rho(df, pred_col="pred_abundance", obs_col="post_abundance"):
    return (df.groupby("recipient")
            .apply(lambda g: spearmanr(g[obs_col], g[pred_col]).statistic,
                   include_groups=False)
            .rename("rho").reset_index())


# ============================================================ load Li data ==
li_all = pd.read_csv(os.path.join(OUT, "AllMetabolicFeatures.csv"))
li = clean(li_all[li_all.timepoint_day == PRIMARY_DAY])

print("=" * 72)
print(f"Li/FAT feature table, day {PRIMARY_DAY}: {li.shape[0]} rows, "
      f"{li.recipient.nunique()} recipients, {li.taxon.nunique()} taxa")
print(li.groupby("fmt_type")[FEATURES].mean().round(3))

# ======================================== (1) Li-native LOO random forest ===
li["loo_pred_fc"] = np.nan
recipients = li.recipient.unique()
imp = {}
for rec in recipients:
    tr = li[li.recipient != rec]
    te = li[li.recipient == rec]
    m = fit_rf(tr)
    li.loc[li.recipient == rec, "loo_pred_fc"] = m.predict(te[FEATURES])
    imp[rec] = m.feature_importances_

li["pred_abundance"] = li.pre_abundance * 2 ** li.loo_pred_fc

rho_li = per_recipient_rho(li).merge(
    li[["recipient", "fmt_type"]].drop_duplicates(), on="recipient")
rho_li["rho_fc"] = [spearmanr(li.loc[li.recipient == r, "log_fold_change"],
                              li.loc[li.recipient == r, "loo_pred_fc"]).statistic
                    for r in rho_li.recipient]

print("\n--- (1) Li-native LOO random forest ---")
print(rho_li.round(3).to_string(index=False))
print(f"\nmedian LOO rho (abundance): {rho_li.rho.median():.3f} +- {rho_li.rho.std():.3f}")
print(f"median LOO rho (fold-change): {rho_li.rho_fc.median():.3f} +- {rho_li.rho_fc.std():.3f}")
print("stitched rho (abundance): "
      f"{spearmanr(li.post_abundance, li.pred_abundance).statistic:.3f}")
print("stitched rho (fold-change): "
      f"{spearmanr(li.log_fold_change, li.loo_pred_fc).statistic:.3f}")

a = rho_li.rho[rho_li.fmt_type == "Allogenic"]
b = rho_li.rho[rho_li.fmt_type == "Autologous"]
print(f"\nAllogenic  rho: {a.median():.3f} +- {a.std():.3f}  (n={len(a)})")
print(f"Autologous rho: {b.median():.3f} +- {b.std():.3f}  (n={len(b)})")
print(f"Wilcoxon p = {mannwhitneyu(a, b).pvalue:.3f}")

imp_df = (pd.DataFrame(imp, index=FEATURES).T
          .agg(["mean", "std"]).T.sort_values("mean", ascending=False))
print("\nPermutation-free impurity importance (mean +- SD across LOO folds):")
print(imp_df.round(4))

# ===================================== (2) IBS -> Li cross-dataset transfer ==
f3 = os.path.join(REPO_DIR, "Figure3")
ibs_dist = pd.read_csv(os.path.join(f3, "processed_data", "nearest_neighbors_t0.csv"))
ibs_dist[["donor", "recipient"]] = ibs_dist.comparison.str.split(" vs ", expand=True)
ibs_div = pd.read_csv(os.path.join(f3, "processed_data", "shannonDiversity.csv")).dropna()
scores = pd.read_csv(os.path.join(f3, "processed_data", "genome_scores.tsv"), sep="\t")
strain = pd.read_csv(os.path.join(f3, "ref", "filtered_strain_taxa.tsv"), sep="\t")
ibs_indep = (scores.merge(strain, left_on="genome_name", right_on="MicrobeID")
             .rename(columns={"genome_score": "met_independence",
                              "NCBI Taxonomy ID": "taxon"})
             .groupby("taxon", as_index=False)["met_independence"].mean())
ibs_abd = pd.read_csv(os.path.join(f3, "processed_data", "pre_post_abundances.csv"))

ibs = (ibs_dist.drop(columns=["sample_id"])
       .merge(ibs_div, on=["donor", "recipient"], how="left")
       .merge(ibs_indep, on="taxon", how="left")
       .merge(ibs_abd, on=["donor", "recipient", "taxon"], how="left"))
ibs["log10_pre_abund"] = np.log10(ibs.pre_abundance)
ibs = clean(ibs)
print("\n" + "=" * 72)
print(f"IBS training table: {ibs.shape[0]} rows, {ibs.recipient.nunique()} recipients")

ibs_model = fit_rf(ibs)
print(f"IBS model in-bag R2: {ibs_model.score(ibs[FEATURES], ibs.log_fold_change):.3f}")

li_x = li.copy()
li_x["xfer_pred_fc"] = ibs_model.predict(li_x[FEATURES])
li_x["xfer_pred_abundance"] = li_x.pre_abundance * 2 ** li_x.xfer_pred_fc

rho_x = per_recipient_rho(li_x, pred_col="xfer_pred_abundance").merge(
    li_x[["recipient", "fmt_type"]].drop_duplicates(), on="recipient")
rho_x["rho_fc"] = [spearmanr(li_x.loc[li_x.recipient == r, "log_fold_change"],
                             li_x.loc[li_x.recipient == r, "xfer_pred_fc"]).statistic
                   for r in rho_x.recipient]

print("\n--- (2) IBS model applied to Li (no Li data in training) ---")
print(rho_x.round(3).to_string(index=False))
print(f"\nmedian transfer rho (abundance): {rho_x.rho.median():.3f} +- {rho_x.rho.std():.3f}")
print(f"median transfer rho (fold-change): {rho_x.rho_fc.median():.3f} +- {rho_x.rho_fc.std():.3f}")
ax = rho_x.rho[rho_x.fmt_type == "Allogenic"]
bx = rho_x.rho[rho_x.fmt_type == "Autologous"]
print(f"Allogenic {ax.median():.3f} vs Autologous {bx.median():.3f}, "
      f"Wilcoxon p = {mannwhitneyu(ax, bx).pvalue:.3f}")

# ------------------------------------------------------------ null models ---
rng = np.random.default_rng(SEED)
null_rho = []
for _ in range(100):
    perm = li.copy()
    perm["pred_abundance"] = perm.pre_abundance * 2 ** rng.permutation(perm.loo_pred_fc.values)
    null_rho.append(per_recipient_rho(perm).rho.median())
null_rho = np.array(null_rho)

# pre-abundance-only baseline: predict no change (log2FC = 0)
li["naive_pred"] = li.pre_abundance
rho_naive = per_recipient_rho(li, pred_col="naive_pred")
print("\n--- nulls / baselines ---")
print(f"shuffled-prediction null rho: {np.median(null_rho):.3f} "
      f"[{np.percentile(null_rho, 2.5):.3f}, {np.percentile(null_rho, 97.5):.3f}]")
print(f"no-change baseline (predict post = pre) rho: {rho_naive.rho.median():.3f}")

# ------------------------------------------------------- timepoint sweep ----
print("\n--- timepoint sweep (Li-native LOO and IBS transfer) ---")
sweep = []
for day in sorted(li_all.timepoint_day.unique()):
    d = clean(li_all[li_all.timepoint_day == day])
    d["loo_pred_fc"] = np.nan
    for rec in d.recipient.unique():
        m = fit_rf(d[d.recipient != rec])
        d.loc[d.recipient == rec, "loo_pred_fc"] = m.predict(d.loc[d.recipient == rec, FEATURES])
    d["pred_abundance"] = d.pre_abundance * 2 ** d.loo_pred_fc
    d["xfer_pred_abundance"] = d.pre_abundance * 2 ** ibs_model.predict(d[FEATURES])
    r_native = per_recipient_rho(d).rho
    r_xfer = per_recipient_rho(d, pred_col="xfer_pred_abundance").rho
    r_naive = per_recipient_rho(d, pred_col="pre_abundance").rho
    sweep.append(dict(day=day, n=len(d),
                      rho_native=r_native.median(), sd_native=r_native.std(),
                      rho_transfer=r_xfer.median(), sd_transfer=r_xfer.std(),
                      rho_nochange=r_naive.median()))
sweep = pd.DataFrame(sweep)
print(sweep.round(3).to_string(index=False))

# ------------------------------------------------------------------ write ---
li.to_csv(os.path.join(OUT, "Li_LOO_predictions_day84.csv"), index=False)
li_x[["recipient", "donor_subject", "fmt_type", "taxon", "pre_abundance",
      "post_abundance", "log_fold_change", "xfer_pred_fc", "xfer_pred_abundance"]].to_csv(
    os.path.join(OUT, "Li_IBStransfer_predictions_day84.csv"), index=False)
rho_li.assign(model="Li-native LOO").to_csv(
    os.path.join(OUT, "Li_native_rho_by_recipient.csv"), index=False)
rho_x.assign(model="IBS transfer").to_csv(
    os.path.join(OUT, "Li_IBStransfer_rho_by_recipient.csv"), index=False)
imp_df.rename_axis("feature").reset_index().to_csv(
    os.path.join(OUT, "Li_feature_importance.csv"), index=False)
sweep.to_csv(os.path.join(OUT, "Li_timepoint_sweep.csv"), index=False)
print(f"\nwrote outputs to {OUT}")
