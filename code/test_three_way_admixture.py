"""
Experiment: Three-way admixture triangle.

Pure EUR, AFR, EAS anchor three corners. Pairwise admixed groups sit on
the triangle edges. A three-way admixed group sits in the interior.
On a PC1 vs PC2 plot this produces a triangle with PC axes rotating to
capture the two dominant admixture axes — demonstrating that PCA
cannot simultaneously resolve 3-way structure without sacrificing
discrimination between populations, and that admixed individuals
falsely appear intermediate rather than as a distinct population.
"""
import runner
from configure import config, PopulationGroup, Population

GROUPS = [
    # Pure anchors
    PopulationGroup(id="pure_EUR",          nsamples=150, populations=[{"EUR": 100}]),
    PopulationGroup(id="pure_AFR",          nsamples=150, populations=[{"AFR": 100}]),
    PopulationGroup(id="pure_EAS",          nsamples=150, populations=[{"EAS": 100}]),
    # Edge midpoints
    PopulationGroup(id="admix_EUR_AFR",     nsamples=100, populations=[{"EUR": 50}, {"AFR": 50}]),
    PopulationGroup(id="admix_EUR_EAS",     nsamples=100, populations=[{"EUR": 50}, {"EAS": 50}]),
    PopulationGroup(id="admix_AFR_EAS",     nsamples=100, populations=[{"AFR": 50}, {"EAS": 50}]),
    # Interior: three-way equal mix
    PopulationGroup(id="admix_three_way",   nsamples=100, populations=[{"EUR": 34}, {"AFR": 33}, {"EAS": 33}]),
]

def main():
    config.mixed_population(GROUPS)
    config.set_phenotype_params(len(GROUPS))
    config.evaluation.metrics.pca = True
    config.evaluation.metrics.gwas = False  # headless Docker segfaults on Pango; use run_gwas_pc_sweep instead
    config.set_output_dir("three_way_admixture")
    config.write()
    # runner.generate_geno(nthreads=3)
    # runner.generate_pheno()
    runner.validate()

if __name__ == "__main__":
    main()
    config.save_config()
