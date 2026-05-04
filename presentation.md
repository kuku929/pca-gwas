---
marp: true
theme: beam
author: Krutarth Patel
title: CS6024 Presentation
paginate: true
---
<!--header: CS6024-->
## **A comparative study of modern population stratification detection methods.**

### *Final Presentation*

*Krutarth Patel*

---

## Background

- **GWAS**(Genome Wide Association Study) tests millions of SNPs for association with a phenotype
- **Population stratification**: systematic ancestry differences lead to false positives in GWAS.
  - Difference in MAF leads to false associations. The associations are caused due to ancestry differences and not due to a direct cause-effect relation.
- **Standard pipline**: include top-*k* PCA eigenvectors as GWAS covariates <sup>[1]</sup>
- **Genomic inflation factor λ**: measures residual stratification

$$\lambda = \frac{\text{median}(\chi^2_{\text{obs}})}{0.456} \quad \text{(ideal: } \lambda = 1\text{)}$$

<small>[1] Price et al. (2006). *Nature Genetics* 38, 904–909.</small>

---

## Is PCA actually robust?

- In practice, *k* is chosen heuristically — k=5 or k=10 by convention
- PCA captures **continuous ancestry gradients**, not discrete clusters <sup>[2]</sup>
- Real biobank cohorts increasingly contain **admixed individuals** who sit between clusters
- PCA output is sensitive to subpopulation sizes, SNP selection, reference panels <sup>[3]</sup>

> Standard evaluations test on clean stratified data.  
> **What happens when assumptions break?**

<small>[2] Novembre & Stephens (2008). *Nature Genetics* 40, 646–649. &nbsp; [3] Elhaik (2022). *Scientific Reports* 12, 10099.</small>

---

## Experimental Design

**Simulator**: HAPNEST — synthetic genotypes from 1KG+HGDP reference panel

| Parameter | Value |
|-----------|-------|
| Chromosome | 2 |
| MAF filter | ≥ 1% |
| Heritability h² | 0.1 uniform *or* per-population (experiment 5) |
| Phenotype model | Single trait, polygenic (0.5% causal SNPs) |

Ground truth is controlled: **any GWAS inflation is purely stratification artifact.**

**GWAS**: plink2 `--glm hide-covar --maf 0.01` with 0 / 2 / 5 / 10 PCs as covariates

---

## Experiments

| # | Name | Design | N |
|---|------|--------|---|
| 1 | Admixture gradient | EUR → 75/25 → 50/50 → 25/75 → AFR | 700 |
| 2 | Three-way admixture | EUR + AFR + EAS anchors + pairwise edges + centre | 850 |
| 3 | Admixed majority overwhelm | 30 EUR + 30 AFR + 500 admixed (EUR/AFR 50/50) | 560 |
| 4 | Admixture sweep | 0% → 25% → 50% → 75% → 100% admixed fraction, fixed N | 700 |
| 5 | Differential heritability | EUR h²=0.3, AFR h²=0.01 — same structure as exp 4 (0%) | 700 |

---

## PCA Structure: Three-Way Admixture

Admixed individuals form **triangle edges** in PC space.

→ [2D Scatter](https://kuku929.github.io/pca-gwas/plots/interactive_pca_scatter_three_way_admixture.html) · [3D (PC1–PC3)](https://kuku929.github.io/pca-gwas/plots/interactive_pca_3d_three_way.html) · [Scree](https://kuku929.github.io/pca-gwas/plots/interactive_scree_three_way_admixture.html)

---

## PCA Structure: Admixed Majority Overwhelm

When 500/560 samples are admixed, **PC axes rotate** to capture variation within the admixed population(to minimize distance).
PCA loses population resolution.

→ [Scatter](https://kuku929.github.io/pca-gwas/plots/interactive_pca_scatter_admixture_overwhelm.html) · [Scree](https://kuku929.github.io/pca-gwas/plots/interactive_scree_admixture_overwhelm.html)

---

## GWAS Inflation Across Experiments

All three cohorts show λ > 1 even with 10 PCs — inflation increases with structural complexity.

**plink2 (PCA correction):**  
→ [Gradient QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_admixture_gradient.html) · [Three-way QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_three_way_admixture.html) · [Overwhelm QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_admixture_overwhelm.html)

**BOLT-LMM (GRM-based):**  
→ [Gradient QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_admixture_gradient_bolt.html) · [Three-way QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_three_way_admixture_bolt.html) · [Overwhelm QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_admixture_overwhelm_bolt.html)

---

## PC Sweep: λ vs Number of PCs

Adding more PCs reduces inflation but **cannot drive λ→1** for admixed cohorts regardless of PC count. Implies that PCA can only account for a few factors. More on that later.  
★ markers show BOLT-LMM reference λ for each experiment.

→ [Open interactive plot](https://kuku929.github.io/pca-gwas/plots/interactive_pc_sweep.html)

---

## Admixture Sweep: λ vs Cohort Composition

Fix N=700 and vary admixed fraction 0→100%.
Inflation peaks at **intermediate admixture** (50%); fully admixed cohort is homogeneous — λ recovers toward 1.  
PCA correction (10 PCs) reduces inflation but it cannot for the population stratification in more admixed cohort.

→ [Open interactive plot](https://kuku929.github.io/pca-gwas/plots/interactive_admixture_sweep.html)

---

## Differential Heritability

**Setup**: pure EUR (h²=0.3) + pure AFR (h²=0.01), N=700 — *identical population structure* to experiment 4.  
PC1 cleanly separates EUR from AFR. Yet λ is high even with 10 PCs.

> PCA corrects mean allele-frequency stratification.
> It cannot correct **variance in effect sizes** across populations.

→ [PC Sweep (with diff heritability)](https://kuku929.github.io/pca-gwas/plots/interactive_pc_sweep.html) · [plink2 QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_diff_heritability.html) · [BOLT-LMM QQ](https://kuku929.github.io/pca-gwas/plots/interactive_qqplot_diff_heritability_bolt.html)

---

## References

1. Price AL, Patterson NJ, Plenge RM, et al. (2006). Principal components analysis corrects for stratification in genome-wide association studies. *Nature Genetics* **38**, 904–909. https://doi.org/10.1038/ng1847

2. Novembre J, Stephens M (2008). Interpreting principal component analyses of spatial population genetic variation. *Nature Genetics* **40**, 646–649. https://doi.org/10.1038/ng.139

3. Elhaik E (2022). Principal Component Analyses (PCA)-based findings in population genetic studies are highly biased and must be reevaluated. *Scientific Reports* **12**, 10099. https://www.nature.com/articles/s41598-022-14395-4

---

## fin.
