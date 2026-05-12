"""
Idea 3: Admixture timing sweep.

Fix EUR/AFR 50/50 admixture proportions; vary when admixture happened (T generations ago).
Hypothesis: λ peaks at intermediate T (~10–20 gens) where long ancestry tracts create LD-driven
inflation that PCA can't resolve. Very recent (T=2) → mosaic genomes. Very ancient (T=200) → panmixia.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import msprime
import numpy as np
from simulate import simulate
from phenotype import stratified_phenotype
from gwas import run_experiment

PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"

TIMING_VALUES = [2, 5, 10, 20, 50, 100, 200]

# EUR=1, ADMIX=0.5, AFR=0 → phenotype correlated with EUR ancestry
POP_SCORES = {"EUR": 1.0, "AFR": 0.0, "ADMIX": 0.5}


def build_demography(T: int) -> tuple:
    """Returns (demography, samples) for admixture at T generations ago."""
    d = msprime.Demography()
    d.add_population(name="EUR",   initial_size=10000)
    d.add_population(name="AFR",   initial_size=14000)
    d.add_population(name="ADMIX", initial_size=10000)
    d.add_population(name="ANC",   initial_size=20000)
    d.add_admixture(time=T, derived="ADMIX", ancestral=["EUR", "AFR"],
                    proportions=[0.5, 0.5])
    d.add_population_split(time=2000, derived=["EUR", "AFR"], ancestral="ANC")
    samples = {"EUR": 200, "AFR": 200, "ADMIX": 300}
    return d, samples


def run(random_seed: int = 42) -> dict:
    """Run timing sweep. Returns {T: {n_pc: lambda}}."""
    all_results = {}
    for T in TIMING_VALUES:
        print(f"\n=== Admixture timing T={T} gen ===")
        out_dir = DATA_DIR / "outputs" / f"msprime_timing_T{T}"
        d, samples = build_demography(T)
        _, ind_names, pop_labels = simulate(
            d, samples, out_dir,
            sequence_length=10e6,
            random_seed=random_seed,
        )
        pheno = stratified_phenotype(pop_labels, POP_SCORES, alpha=0.8, random_seed=random_seed)
        lambdas = run_experiment(out_dir, ind_names, pop_labels, pheno)
        all_results[T] = lambdas

    summary_path = DATA_DIR / "outputs" / "msprime_timing_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n=> Summary: {summary_path}")
    _print_table(all_results)
    return all_results


def _print_table(results: dict):
    pc_counts = sorted({k for v in results.values() for k in v})
    header = f"{'T':>5}" + "".join(f"  PC{n:>2}" for n in pc_counts)
    print(header)
    for T in sorted(results):
        row = f"{T:>5}" + "".join(
            f"  {results[T].get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        print(row)


if __name__ == "__main__":
    run()
