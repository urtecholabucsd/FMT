# FMT metabolic modeling analyses

This repository contains analysis scripts and processed data supporting the current [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.05.15.725570v1), *Metagenome-scale Modeling to Assess Microbiome Metabolic Complementarity for Precision Microbiota Transplantation Therapies*. The manuscript link will be updated after publication.

## Getting started

Clone the repository and create the analysis environment with Conda:

```bash
git clone https://github.com/urtecholabucsd/FMT.git
cd FMT
conda env create -f environment.yml
conda activate fmt-analysis
```

The R analyses use a separate R installation. Required R packages are declared with `library()` near the beginning of each document and are not currently captured in `environment.yml`.

If you are interested in building the model databases, please refer to the installation guide for [q2-micom plugin](https://micom-dev.github.io/q2-micom/).

## Repository organization

- `Figure1/`–`Figure4/`: processed data and scripts used for the main figures.
- `goll_et.al_2020_IBS_FMT/`: preprocessing and modeling of the human IBS FMT cohort.
- `li_et.al_2016_FMT/`: preprocessing and processed results for the Li et al. 2016 cohort.
- `urtecho_et.al_2024_mouse_FMT/`: preprocessing and modeling of the mouse FMT study.
- `OpenBiome_Predictions/`: raw inputs, processed predictions, and scripts for the OpenBiome donor analysis.
- `utils/`: shared metadata, plotting resources, and helper functions.

Large files are not committed to the repository. They can be downloaded from the corresponding data repositories or regenerated using the commands documented below.

## Source datasets

Raw sequencing data are not duplicated in this repository. The analyses use data from the following published studies and public archives:

| Analysis directory | Source study | Raw-data accession |
| --- | --- | --- |
| `goll_et.al_2020_IBS_FMT/` | Goll et al. (2020), [*Effects of fecal microbiota transplantation in subjects with irritable bowel syndrome are mirrored by changes in gut microbiome*](https://www.tandfonline.com/doi/full/10.1080/19490976.2020.1794263) | [ENA PRJEB36140](https://www.ebi.ac.uk/ena/browser/view/PRJEB36140) |
| `li_et.al_2016_FMT/` | Li et al. (2016), [*Durable coexistence of donor and recipient strains after fecal microbiota transplantation*](https://www.science.org/doi/10.1126/science.aad8852) | [ENA/NCBI BioProject PRJEB12357](https://www.ncbi.nlm.nih.gov/bioproject/PRJEB12357) |
| `urtecho_et.al_2024_mouse_FMT/` | Urtecho et al. (2024), [*Spatiotemporal dynamics during niche remodeling by super-colonizing microbiota in the mammalian gut*](https://www.cell.com/cell-systems/fulltext/S2405-4712(24)00304-1) | [NCBI BioProject PRJNA1028308](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA1028308) |
| `OpenBiome_Predictions/` | Poyet et al. (2019), [*A library of human gut bacterial isolates paired with longitudinal multiomics data enables mechanistic microbiome research*](https://doi.org/10.1038/s41591-019-0559-3) | [NCBI BioProject PRJNA544527](https://www.ncbi.nlm.nih.gov/bioproject/544527) |

## Software and model resources

The study uses the following external software, models, and media repositories:

- [MICOM](https://github.com/micom-dev/micom) for community metabolic modeling and flux balance analysis.
- [CarveMe](https://github.com/cdanielmachado/carveme) for mouse genome-scale metabolic model reconstruction from metagenome-assembled genomes.
- [AGORA](https://github.com/VirtualMetabolicHuman/AGORA) for computing metabolic independence from metagenome-assembled genomes.
- [MICOM media](https://github.com/micom-dev/media) for growth media definitions used in MICOM.
- [q2-micom](https://micom-dev.github.io/q2-micom/) for building model database at user-defined taxonomic levels for MICOM.

## Citing this work

If you use this repository, please cite the current preprint:

> Zhang Z, Holton M, Ferrer DM, Tripp AD, Richter A, Dixit PD, Urtecho G. *Metagenome-scale Modeling to Assess Microbiome Metabolic Complementarity for Precision Microbiota Transplantation Therapies*. bioRxiv (2026). [https://doi.org/10.64898/2026.05.15.725570](https://www.biorxiv.org/content/10.64898/2026.05.15.725570v1)

Please also cite MICOM, CarveMe, AGORA, and the source dataset publications relevant to the analysis being reused.

## License

The original code in this repository is available under the [MIT License](LICENSE). Source datasets and third-party software remain subject to their own licenses and terms of use.
