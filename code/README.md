# PCA Robustness under Admixture — GWAS Stratification Study

Testing failure modes of PCA-based population stratification correction in GWAS using HAPNEST-simulated genotypes.

## Overview

Principal Component Analysis (PCA) is the standard method for correcting population stratification in genome-wide association studies (GWAS). This project stress-tests PCA against a range of admixed cohort designs to expose its failure modes, quantified via the genomic inflation factor **λ**.

All experiments use synthetic genotypes from the HAPNEST simulator (1KG+HGDP reference, chromosome 2). Ground truth heritability is fixed so that any GWAS inflation is a pure stratification artifact. BOLT-LMM (GRM-based correction) is run alongside as a comparison method.

## Experiments

| # | Name | Design | N |
|---|------|--------|---|
| 1 | Admixture gradient | EUR → 75/25 → 50/50 → 25/75 → AFR | 700 |
| 2 | Three-way admixture | EUR + AFR + EAS anchors + pairwise edges + centre | 850 |
| 3 | Admixed majority overwhelm | 30 EUR + 30 AFR + 500 admixed (EUR/AFR 50/50) | 560 |
| 4 | Admixture sweep | 0% → 25% → 50% → 75% → 100% admixed fraction, fixed N | 700 |
| 5 | Differential heritability | EUR h²=0.3, AFR h²=0.01 — same structure as exp 4 (0%) | 700 |

## Requirements

- Docker container named `hapnest` (running HAPNEST)
- Python ≥ 3.13 with `uv`
- BOLT-LMM v2.5 binary (see below)

```bash
uv sync
```

### BOLT-LMM installation

```bash
cd ..   # project root
wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz
tar xzf BOLT-LMM_v2.5.tar.gz
rm BOLT-LMM_v2.5.tar.gz
```

### Docker setup (one-time)

```bash
docker exec hapnest chmod 755 /root
```

## Data

Input reference data goes in `../data/inputs/` (not tracked in git — download separately via HAPNEST docs).
Generated outputs go in `../data/outputs/test_<name>/`.

## Running

```bash
# Generate all experiments (each ~5-10 min)
uv run python test_admixture_gradient.py
uv run python test_three_way_admixture.py
uv run python test_admixture_overwhelm.py
uv run python test_admixture_sweep.py
uv run python test_diff_heritability.py

# BOLT-LMM on all experiments (~5 min each)
uv run python test_bolt_lmm.py

# Visualisations (fast)
uv run python ../viz_pca.py admixture_gradient three_way_admixture admixture_overwhelm diff_heritability
uv run python ../viz_gwas.py
uv run python ../viz_gwas.py --bolt admixture_gradient three_way_admixture admixture_overwhelm diff_heritability
uv run python ../viz_pc_sweep.py
```

## Key files

| File | Purpose |
|------|---------|
| `configure.py` | Pydantic config model + helpers (`mixed_population`, `set_phenotype_params`) |
| `runner.py` | HAPNEST Docker wrapper + plink2 PC sweep + GWAS λ computation |
| `test_*.py` | Per-experiment scripts |
| `../viz_pca.py` | PCA scatter, scree, 3D plots |
| `../viz_gwas.py` | QQ plots (plink2 + BOLT-LMM) |
| `../viz_pc_sweep.py` | λ vs n_PCs sweep + admixture sweep charts |
| `../presentation.md` | Marp presentation |

## References

1. Price et al. (2006). *Nature Genetics* 38, 904–909. https://doi.org/10.1038/ng1847
2. Novembre & Stephens (2008). *Nature Genetics* 40, 646–649. https://doi.org/10.1038/ng.139
3. Elhaik (2022). *Scientific Reports* 12, 10099. https://www.nature.com/articles/s41598-022-14395-4
