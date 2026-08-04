#!/usr/bin/env python3
"""05_figures.py - summary panels for the Li/FAT replication."""

import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BASE = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(BASE, "data", "supp_5")
FIGS = os.path.join(BASE, "figures")
os.makedirs(FIGS, exist_ok=True)
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9,
                     "axes.spines.top": False, "axes.spines.right": False})

li = pd.read_csv(os.path.join(OUT, "Li_LOO_predictions_day84.csv"))
xf = pd.read_csv(os.path.join(OUT, "Li_IBStransfer_predictions_day84.csv"))
abl = pd.read_csv(os.path.join(OUT, "feature_ablation.csv"))
grp = pd.read_csv(os.path.join(OUT, "Li_by_group_and_timepoint.csv"))
base = pd.read_csv(os.path.join(OUT, "baseline_control_IBS_vs_Li.csv"))

fig, ax = plt.subplots(2, 2, figsize=(9, 7))

# A. observed vs predicted post-FMT abundance, Li-native LOO, day 84
a = ax[0, 0]
c = {"Allogenic": "#2c7fb8", "Autologous": "#d95f02"}
for g, d in li.groupby("fmt_type"):
    a.scatter(np.log10(d.post_abundance + 1e-6), np.log10(d.pred_abundance + 1e-6),
              s=6, alpha=.35, color=c[g], label=g, edgecolors="none")
lim = [1, 6]
a.plot(lim, lim, "--", color="grey", lw=.8)
a.set(xlabel="Observed log10(post-FMT abundance)",
      ylabel="Predicted log10(post-FMT abundance)",
      title="A  Li-native LOO, day 84 (rho = 0.78)")
a.legend(frameon=False, fontsize=8)

# B. fold-change rho: this is the metric that actually tests the model
b = ax[0, 1]
labels = ["IBS\nmodel", "IBS\nno-change", "Li\nmodel", "Li\nno-change"]
vals = [base.loc[(base.dataset.str.startswith("IBS")) &
                 (base.metric == "rho fold-change, model"), "median"].item(),
        0.0,
        base.loc[(base.dataset.str.startswith("Li")) &
                 (base.metric == "rho fold-change, model"), "median"].item(),
        0.0]
abd = [base.loc[(base.dataset.str.startswith("IBS")) &
                (base.metric == "rho abundance, model"), "median"].item(),
       base.loc[(base.dataset.str.startswith("IBS")) &
                (base.metric == "rho abundance, no-change baseline"), "median"].item(),
       base.loc[(base.dataset.str.startswith("Li")) &
                (base.metric == "rho abundance, model"), "median"].item(),
       base.loc[(base.dataset.str.startswith("Li")) &
                (base.metric == "rho abundance, no-change baseline"), "median"].item()]
x = np.arange(4)
b.bar(x - .2, abd, .38, color="#bdbdbd", edgecolor="k", lw=.6, label="abundance-level rho")
b.bar(x + .2, vals, .38, color="#2c7fb8", edgecolor="k", lw=.6, label="fold-change rho")
b.set_xticks(x); b.set_xticklabels(labels, fontsize=8)
b.axhline(0, color="k", lw=.6)
b.set(ylabel="Spearman rho", ylim=(-0.05, 1.0),
      title="B  Abundance-level rho is not model-informative")
b.legend(frameon=False, fontsize=8, loc="upper right")

# C. fold-change rho across timepoints and fmt_type
cax = ax[1, 0]
for g, m in [("Allogenic", "o"), ("Autologous", "s"), ("all", "^")]:
    d = grp[grp.group == g].sort_values("day")
    cax.errorbar(d.day, d.rho_fc, yerr=d.sd_fc, marker=m, capsize=3, lw=1.2,
                 label=g, color=c.get(g, "grey"))
cax.axhline(0, color="k", lw=.6)
cax.axhline(0.535, ls="--", color="firebrick", lw=1)
cax.text(45, 0.55, "IBS LOO (6 months)", color="firebrick", fontsize=8)
cax.set(xlabel="Days post-FMT", ylabel="LOO fold-change Spearman rho",
        title="C  Predictive signal decays after transplant")
cax.legend(frameon=False, fontsize=8)

# D. feature ablation
d = ax[1, 1]
order = ["full (4 features)", "log10_pre_abund only", "metabolic only (no pre-abund)"]
w = .38
for i, (ds, col) in enumerate([("IBS", "#2c7fb8"), ("Li day42", "#d95f02")]):
    v = [abl.loc[(abl.dataset == ds) & (abl.features == f), "rho_fc"].item() for f in order]
    e = [abl.loc[(abl.dataset == ds) & (abl.features == f), "sd_fc"].item() for f in order]
    d.bar(np.arange(3) + (i - .5) * w, v, w, yerr=e, capsize=3,
          color=col, edgecolor="k", lw=.6, label=ds)
d.set_xticks(range(3))
d.set_xticklabels(["Full\n(4 features)", "Pre-abundance\nonly", "Metabolic\nonly"], fontsize=8)
d.axhline(0, color="k", lw=.6)
d.set(ylabel="LOO fold-change Spearman rho", title="D  Feature ablation")
d.legend(frameon=False, fontsize=8)

fig.tight_layout()
fig.savefig(os.path.join(FIGS, "Li_replication_summary.png"), dpi=300, facecolor="white")
fig.savefig(os.path.join(FIGS, "Li_replication_summary.pdf"), facecolor="white")
print("wrote", os.path.join(FIGS, "Li_replication_summary.png"))
