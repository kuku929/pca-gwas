# PCA Robustness under Admixture — GWAS Stratification Study

Testing failure modes of PCA-based population stratification correction in GWAS using HAPNEST-simulated genotypes.

## Overview

Principal Component Analysis (PCA) is the standard method for correcting population stratification in genome-wide association studies (GWAS). This project stress-tests PCA against a range of admixed cohort designs to expose its failure modes, quantified via the genomic inflation factor **λ**.

All experiments use synthetic genotypes from the HAPNEST simulator (1KG+HGDP reference, chromosome 2). Ground truth heritability is fixed so that any GWAS inflation is a pure stratification artifact. BOLT-LMM (GRM-based correction) is run alongside as a comparison method.

## Interactive Plots

All plots are hosted at [kuku929.github.io/pca-gwas](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main).

| Plot | Link |
|------|------|
| PC Sweep: λ vs n_PCs | [interactive_pc_sweep](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main/plots/interactive_pc_sweep.html) |
| Admixture Sweep: λ vs cohort composition | [interactive_admixture_sweep](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main/plots/interactive_admixture_sweep.html) |
| Three-way PCA scatter (2D) | [pca_scatter](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main/plots/interactive_pca_scatter_three_way_admixture.html) |
| Three-way PCA scatter (3D) | [pca_3d](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main/plots/interactive_pca_3d_three_way.html) |
| Admixed majority overwhelm PCA | [pca_overwhelm](https://htmlpreview.github.io/?https://raw.githubusercontent.com/kuku929/pca-gwas/main/plots/interactive_pca_scatter_admixture_overwhelm.html) |

## Experiments

| # | Name | Design | N |
|---|------|--------|---|
| 1 | Admixture gradient | EUR → 75/25 → 50/50 → 25/75 → AFR | 700 |
| 2 | Three-way admixture | EUR + AFR + EAS anchors + pairwise edges + centre | 850 |
| 3 | Admixed majority overwhelm | 30 EUR + 30 AFR + 500 admixed (EUR/AFR 50/50) | 560 |
| 4 | Admixture sweep | 0% → 25% → 50% → 75% → 100% admixed fraction, fixed N | 700 |
| 5 | Differential heritability | EUR h²=0.3, AFR h²=0.01 — same structure as exp 4 (0%) | 700 |

## Code structure

```
code/           Python experiment scripts and HAPNEST runner
viz_pca.py      PCA scatter, scree, 3D plots
viz_gwas.py     QQ plots (plink2 + BOLT-LMM)
viz_pc_sweep.py λ vs n_PCs + admixture sweep charts
presentation.md Marp presentation
plots/          Pre-generated interactive HTML plots
```

See `code/README.md` for full setup and run instructions.

## References

1. Price et al. (2006). *Nature Genetics* 38, 904–909. https://doi.org/10.1038/ng1847
2. Novembre & Stephens (2008). *Nature Genetics* 40, 646–649. https://doi.org/10.1038/ng.139
3. Elhaik (2022). *Scientific Reports* 12, 10099. https://www.nature.com/articles/s41598-022-14395-4
