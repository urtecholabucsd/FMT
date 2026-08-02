# FMT metabolic modeling analyses

This repository contains the analysis code and processed data supporting the current [bioRxiv preprint](https://www.biorxiv.org/content/10.64898/2026.05.15.725570v1). The manuscript link and citation will be updated after publication.

## Repository organization

- `Figure1/`–`Figure4/`: processed data and scripts used for the main figures.
- `goll_et.al_2020_IBS_FMT/`: preprocessing and modeling of the human IBS FMT cohort.
- `li_et.al_2016_FMT/`: preprocessing and processed results for the Li et al. 2016 cohort.
- `urtecho_et.al_2024_mouse_FMT/`: preprocessing and modeling of the mouse FMT study.
- `OpenBiome_Predictions/`: raw inputs, processed predictions, and scripts for the OpenBiome donor analysis.
- `utils/`: shared metadata and plotting resources.

Analysis documents use paths relative to their locations. No formal package installation procedure is provided for this repository. Large generated model archives are kept outside version control and can be recreated using the commands that will be documented below.

## Software and model resources

The study uses the following external software, models, and media repositories:

- [MICOM](https://github.com/micom-dev/micom) for community metabolic modeling.
- [CarveMe](https://github.com/cdanielmachado/carveme) for genome-scale metabolic model reconstruction.
- [AGORA](https://github.com/VirtualMetabolicHuman/AGORA) as the microbial genome-scale model resource.
- [AGORA version 1.03](https://github.com/VirtualMetabolicHuman/AGORA/tree/master/CurrentVersion/AGORA_1_03) for the model collection used in this study.
- [MICOM media](https://github.com/micom-dev/media) for diet and growth media definitions.

## Reproduction commands

Exact commands used to build models, run simulations, and reproduce the analyses will be added here.

### Model construction

```text
TODO: add exact model-construction commands.
```

### Community simulations

```text
TODO: add exact MICOM simulation commands.
```

### Figure generation

```text
TODO: add exact notebook and R Markdown rendering commands.
```
