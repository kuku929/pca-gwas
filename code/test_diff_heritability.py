"""
Controlled experiment: differential heritability across populations.

Same population structure as the 0% admixture sweep baseline (pure EUR + pure AFR,
N=700) — but heritability differs: EUR h²=0.3, AFR h²=0.01.

PCA with 1 PC perfectly separates EUR from AFR. Yet GWAS lambda stays elevated
because PCA corrects mean allele-frequency stratification, not variance in effect
sizes across populations. Any residual inflation after 10 PCs is attributable
purely to genetic architecture heterogeneity, not ancestry confounding.

Compare against admixture_sweep_0pct (same structure, uniform h²=0.1) to isolate
the heritability effect.
"""
import runner
from configure import config, PopulationGroup

GROUPS = [
    PopulationGroup(id="pure_EUR", nsamples=350, populations=[{"EUR": 100}]),
    PopulationGroup(id="pure_AFR", nsamples=350, populations=[{"AFR": 100}]),
]

def main():
    config.mixed_population(GROUPS)
    config.set_phenotype_params(len(GROUPS), heritability=[0.3, 0.01])
    config.evaluation.metrics.pca  = True
    config.evaluation.metrics.gwas = False
    config.set_output_dir("diff_heritability")
    config.write()

    runner.generate_geno(nthreads=3)
    runner.generate_pheno()
    runner.validate()

if __name__ == "__main__":
    main()
    config.save_config()
