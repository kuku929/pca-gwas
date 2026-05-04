"""
Experiment: EUR-AFR admixture gradient.

5 groups with increasing AFR ancestry (0% → 25% → 50% → 75% → 100%).
Admixed individuals form a gradient on PC1 between the two pure populations,
demonstrating that PCA axes capture ancestry proportion, not just
discrete population membership — confounding GWAS stratification correction.
"""
import runner
from configure import config, PopulationGroup, Population

GROUPS = [
    PopulationGroup(id="pure_EUR",     nsamples=200, populations=[{"EUR": 100}]),
    PopulationGroup(id="admix_75EUR",  nsamples=100, populations=[{"EUR": 75}, {"AFR": 25}]),
    PopulationGroup(id="admix_50EUR",  nsamples=100, populations=[{"EUR": 50}, {"AFR": 50}]),
    PopulationGroup(id="admix_25EUR",  nsamples=100, populations=[{"EUR": 25}, {"AFR": 75}]),
    PopulationGroup(id="pure_AFR",     nsamples=200, populations=[{"AFR": 100}]),
]

def main():
    config.mixed_population(GROUPS)
    config.set_phenotype_params(len(GROUPS))
    config.evaluation.metrics.pca = True
    config.evaluation.metrics.gwas = False # headless Docker segfaults on Pango; use run_gwas_pc_sweep instead
    config.set_output_dir("admixture_gradient")
    config.write()
    runner.generate_geno(nthreads=3)
    runner.generate_pheno()
    runner.validate()

if __name__ == "__main__":
    main()
    config.save_config()
