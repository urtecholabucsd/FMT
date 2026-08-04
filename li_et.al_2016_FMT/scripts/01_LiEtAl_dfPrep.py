#!/usr/bin/env python3
"""
01_LiEtAl_dfPrep.py

Prepare model-input dataframes for the Li et al. FAT cohort (PRJEB46778),
mirroring the IBS preprocessing workflow used for Figure 3.

Design decisions (confirmed with GU):
  - post-FMT timepoint = day 42 (primary); day 2/14/84 also written for the sweep
  - all 10 recipients retained (5 Allogenic, 5 Autologous); fmt_type replaces
    clinical_response as the stratifying variable
  - donor FAT_DON_11 has 3 sequencing replicates; normalized abundances and
    metabolic distances are averaged across them

Writes the processed inputs for Supplementary Figure 5 to
Figure2/data/supp_5/.
"""

import os
import numpy as np
import pandas as pd

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
COHORT_DIR = os.path.dirname(SCRIPT_DIR)
REPO_DIR = os.path.dirname(COHORT_DIR)
RAW_DIR = os.path.join(COHORT_DIR, "raw")
PROCESSED_DIR = os.path.join(COHORT_DIR, "processed")
OUT = os.path.join(REPO_DIR, "Figure2", "data", "supp_5")
os.makedirs(OUT, exist_ok=True)

POST_DAYS = [2, 14, 42, 84]

# ---------------------------------------------------------------- metadata --
meta = pd.read_excel(os.path.join(RAW_DIR, "PRJEB46778_categorized_samples_ALL.xlsx"))
pairs = pd.read_excel(os.path.join(RAW_DIR, "FAT_recipient_donor_pairs.xlsx"))

counts = pd.read_csv(os.path.join(RAW_DIR, "LiEtAl_agora_mapped_counts.csv"))
samples = [c for c in counts.columns if c != "taxon_id"]

meta = meta[meta.run_accession.isin(samples)].copy()

# recipient pre / post run accessions
rec_pre = (meta[meta.sample_role == "Recipient_baseline_preFMT"]
           .set_index("subject")["run_accession"].to_dict())
rec_post = {d: (meta[(meta.sample_role == "Recipient_postFMT") & (meta.timepoint_day == d)]
                .set_index("subject")["run_accession"].to_dict())
            for d in POST_DAYS}

# donor subject -> list of replicate runs
donor_runs = (meta[meta.sample_role == "Donor"]
              .groupby("subject")["run_accession"].apply(list).to_dict())
# autologous donors are the subject's own baseline sample
for _, r in pairs.iterrows():
    if r.fmt_type == "Autologous":
        donor_runs.setdefault(r.donor, [rec_pre[r.recipient]])

mapping = pairs.rename(columns={"recipient": "recipient", "donor": "donor_subject"}).copy()
mapping["recipient_pre"] = mapping.recipient.map(rec_pre)
for d in POST_DAYS:
    mapping[f"recipient_post_d{d}"] = mapping.recipient.map(rec_post[d])
mapping["donor_replicates"] = mapping.donor_subject.map(lambda s: ";".join(donor_runs[s]))
mapping.to_csv(os.path.join(OUT, "DonorRecipientMapping.csv"), index=False)

# ------------------------------------------------- normalize count matrix ---
mat = counts.set_index("taxon_id").astype(float) + 1.0   # pseudocount, as in IBS prep
col_sums = mat.sum(axis=0)
mat_norm = mat.divide(col_sums, axis=1) * np.median(col_sums)

# ------------------------------------------------------ Shannon diversity ---
def shannon(v):
    p = v / v.sum()
    p = p[p > 0]
    return float(-(p * np.log(p)).sum())

shannon_run = pd.Series({s: shannon(mat.loc[:, s].values) for s in samples}, name="shannonD")
shannon_run.rename_axis("sample").reset_index().to_csv(
    os.path.join(OUT, "ShannonDF_All.csv"), index=False)

# donor-level profiles: mean normalized abundance across replicates
donor_profile = {sub: mat_norm[runs].mean(axis=1) for sub, runs in donor_runs.items()}
donor_shannon = {sub: shannon(p.values) for sub, p in donor_profile.items()}

div = mapping[["donor_subject", "recipient", "fmt_type"]].copy()
div["donor_shannonD"] = div.donor_subject.map(donor_shannon)
div["recipient_shannonD"] = mapping.recipient_pre.map(shannon_run)
div["diversity_ratio"] = div.donor_shannonD / div.recipient_shannonD
div.to_csv(os.path.join(OUT, "shannonDiversity.csv"), index=False)

# ------------------------------------------------------------- abundances ---
long = (mat_norm.reset_index()
        .melt(id_vars="taxon_id", var_name="sample", value_name="normalized_count"))
long.to_csv(os.path.join(OUT, "TaxonCounts_all.csv"), index=False)

donor_abd = (pd.DataFrame(donor_profile)
             .rename_axis("taxon").reset_index()
             .melt(id_vars="taxon", var_name="donor_subject", value_name="donor_abundance"))
donor_abd.to_csv(os.path.join(OUT, "donor_abundances.csv"), index=False)

pre_post = []
for _, r in mapping.iterrows():
    pre = mat_norm[r.recipient_pre]
    for d in POST_DAYS:
        post = mat_norm[r[f"recipient_post_d{d}"]]
        pre_post.append(pd.DataFrame({
            "recipient": r.recipient,
            "donor_subject": r.donor_subject,
            "fmt_type": r.fmt_type,
            "timepoint_day": d,
            "taxon": mat_norm.index,
            "pre_abundance": pre.values,
            "post_abundance": post.values,
            "log_fold_change": np.log2(post.values / pre.values),
        }))
pre_post = pd.concat(pre_post, ignore_index=True)
pre_post.to_csv(os.path.join(OUT, "pre_post_abundances.csv"), index=False)

# ------------------------------------------- metabolic interaction distance --
nn = pd.read_csv(os.path.join(PROCESSED_DIR, "nearest_neighbors.csv"))
nn[["s1", "s2"]] = nn.comparison.str.split(" vs ", expand=True)
# comparisons are stored once per unordered pair; index both orientations
nn_key = pd.concat([
    nn.assign(a=nn.s1, b=nn.s2),
    nn.assign(a=nn.s2, b=nn.s1),
], ignore_index=True).drop_duplicates(subset=["a", "b", "taxon"])

dist_rows = []
for _, r in mapping.iterrows():
    reps = donor_runs[r.donor_subject]
    sub = nn_key[(nn_key.a.isin(reps)) & (nn_key.b == r.recipient_pre)]
    if sub.empty:
        raise RuntimeError(f"no distance found for {r.donor_subject} vs {r.recipient}")
    agg = (sub.groupby("taxon")["average_distance"].mean().reset_index())
    agg["recipient"] = r.recipient
    agg["donor_subject"] = r.donor_subject
    agg["n_donor_replicates_used"] = sub.a.nunique()
    dist_rows.append(agg)
dist = pd.concat(dist_rows, ignore_index=True)
dist.to_csv(os.path.join(OUT, "metabolic_distance_t0.csv"), index=False)

# ------------------------------------------------- metabolic independence ---
scores = pd.read_csv(os.path.join(REPO_DIR, "Figure3", "processed_data",
                                  "genome_scores.tsv"), sep="\t")
strain = pd.read_csv(os.path.join(REPO_DIR, "Figure3", "ref",
                                  "filtered_strain_taxa.tsv"), sep="\t")
indep = (scores.merge(strain, left_on="genome_name", right_on="MicrobeID")
         .rename(columns={"genome_score": "met_independence",
                          "NCBI Taxonomy ID": "taxon"})[["taxon", "met_independence"]])
indep = indep.groupby("taxon", as_index=False)["met_independence"].mean()
indep.to_csv(os.path.join(OUT, "metabolic_independence.csv"), index=False)

# --------------------------------------------------- assemble feature table --
feat = (dist
        .merge(div[["donor_subject", "recipient", "fmt_type", "diversity_ratio",
                    "donor_shannonD", "recipient_shannonD"]],
               on=["donor_subject", "recipient"], how="left")
        .merge(indep, on="taxon", how="left")
        .merge(pre_post, on=["recipient", "donor_subject", "taxon"], how="left")
        .merge(donor_abd, on=["donor_subject", "taxon"], how="left"))
feat = feat.drop(columns=["fmt_type_y"]).rename(columns={"fmt_type_x": "fmt_type"})
feat["log10_pre_abund"] = np.log10(feat.pre_abundance)
feat.to_csv(os.path.join(OUT, "AllMetabolicFeatures.csv"), index=False)

print(f"mapping: {len(mapping)} donor-recipient pairs")
print(f"feature table: {feat.shape[0]} rows, {feat.taxon.nunique()} taxa, "
      f"{feat.recipient.nunique()} recipients, {len(POST_DAYS)} timepoints")
print(feat.groupby(['fmt_type', 'timepoint_day']).size())
print("missing met_independence:", feat.met_independence.isna().sum())
