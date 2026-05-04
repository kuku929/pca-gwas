"""
Cohort composition sweep: 5 separate GWAS runs with increasing proportions
of admixed (50/50 EUR-AFR) individuals in the cohort (0 → 25 → 50 → 75 → 100%).
Total N fixed at 700 per run.

Key question: does PCA correction (10 PCs) degrade monotonically as the cohort
becomes more admixed? Expected: lambda rises then drops at 100% (fully admixed
cohort is homogeneous — no stratification to correct for). Peak inflation at
intermediate admixture proportions where admixed individuals break PC structure.

Each run is independent and takes ~5-10 minutes.
"""
import runner
from configure import config, PopulationGroup

SWEEP = [
    ("0pct", [
        PopulationGroup(id="pure_EUR", nsamples=350, populations=[{"EUR": 100}]),
        PopulationGroup(id="pure_AFR", nsamples=350, populations=[{"AFR": 100}]),
    ]),
    ("25pct", [
        PopulationGroup(id="pure_EUR", nsamples=262, populations=[{"EUR": 100}]),
        PopulationGroup(id="pure_AFR", nsamples=263, populations=[{"AFR": 100}]),
        PopulationGroup(id="admixed",  nsamples=175, populations=[{"EUR": 50}, {"AFR": 50}]),
    ]),
    ("50pct", [
        PopulationGroup(id="pure_EUR", nsamples=175, populations=[{"EUR": 100}]),
        PopulationGroup(id="pure_AFR", nsamples=175, populations=[{"AFR": 100}]),
        PopulationGroup(id="admixed",  nsamples=350, populations=[{"EUR": 50}, {"AFR": 50}]),
    ]),
    ("75pct", [
        PopulationGroup(id="pure_EUR", nsamples=87,  populations=[{"EUR": 100}]),
        PopulationGroup(id="pure_AFR", nsamples=88,  populations=[{"AFR": 100}]),
        PopulationGroup(id="admixed",  nsamples=525, populations=[{"EUR": 50}, {"AFR": 50}]),
    ]),
    ("100pct", [
        PopulationGroup(id="admixed",  nsamples=700, populations=[{"EUR": 50}, {"AFR": 50}]),
    ]),
]


def run_sweep_point(label: str, groups: list):
    suffix = f"admixture_sweep_{label}"
    config.mixed_population(groups)
    config.set_phenotype_params(len(groups))
    config.evaluation.metrics.pca  = True
    config.evaluation.metrics.gwas = False  # headless Docker segfaults on Pango
    config.set_output_dir(suffix)
    config.write()

    runner.generate_geno(nthreads=3)
    runner.generate_pheno()
    runner.validate()
    runner.run_gwas_pc_sweep(
        output_dir=config.filepaths.general.output_dir,
        prefix_template=config.filepaths.general.output_prefix,
        chromosome=config.global_parameters.chromosome,
        pc_counts=[0, 10],   # no-correction vs 10-PC correction
    )
    config.save_config()


def main():
    for label, groups in SWEEP:
        print(f"\n{'='*55}")
        print(f"=== SWEEP {label}  ({len(groups)} groups, {sum(g.nsamples for g in groups)} samples) ===")
        print(f"{'='*55}")
        run_sweep_point(label, groups)


if __name__ == "__main__":
    main()
