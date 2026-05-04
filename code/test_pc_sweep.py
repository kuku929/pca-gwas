"""
PC count sweep: run plink2 GWAS with 0, 2, 5, 10 PCs as covariates on
existing gradient experiment data. No data regeneration needed.

Shows that increasing PC count does NOT drive lambda to 1 for admixed cohorts —
PCA correction is fundamentally limited regardless of how many PCs are used.
"""
import runner

def main():
    runner.run_gwas_pc_sweep(
        output_dir="data/outputs/test_admixture_gradient",
        prefix_template="test_chr-{chromosome}",
        chromosome=2,
        pc_counts=[0, 2, 5, 10],
    )

if __name__ == "__main__":
    main()
