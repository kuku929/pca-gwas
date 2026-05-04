"""
Experiment: Admixed majority overwhelms PCA.

Two small pure-population groups (30 each) surrounded by a large admixed
cohort (500). When admixed samples dominate the variance, PC1 and PC2
rotate to capture within-admixed variation — the pure-population clusters
collapse toward the center and lose separation. This shows PCA-based
stratification correction fails when cohort composition is unbalanced,
a realistic scenario in biobank studies with many admixed participants.
"""
import runner
from configure import config, PopulationGroup, Population

GROUPS = [
    PopulationGroup(id="pure_EUR",         nsamples=30,  populations=[{"EUR": 100}]),
    PopulationGroup(id="pure_AFR",         nsamples=30,  populations=[{"AFR": 100}]),
    PopulationGroup(id="admixed_majority", nsamples=500, populations=[{"EUR": 50}, {"AFR": 50}]),
]

def main():
    config.mixed_population(GROUPS)
    config.set_phenotype_params(len(GROUPS))
    config.evaluation.metrics.pca = True
    config.evaluation.metrics.gwas = False  # headless Docker segfaults on Pango; use run_gwas_pc_sweep instead
    config.set_output_dir("admixture_overwhelm")
    config.write()
    # runner.generate_geno(nthreads=3)
    runner.generate_pheno()
    runner.validate()

if __name__ == "__main__":
    main()
    config.save_config()
