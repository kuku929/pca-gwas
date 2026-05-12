"""
Idea 1: Isolation by distance (1D stepping-stone model).

K populations arranged on a line; migration only between neighbours.
Phenotype = geographic position (pop_index / (K-1)).
PCA eigenvectors on clinal variation are sinusoidal (Novembre et al. 2008) —
not useful as stratification covariates. Even 20 PCs cannot fully correct λ
because the confounding is infinite-dimensional (continuous geography).

Also tests: what happens if we sample only from the two endpoints (typical
biobank case — two discrete populations) vs sampling the full gradient.
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

K = 16          # number of populations along the stepping stone
N_EACH = 44     # individuals per population → total ~700
M = 1e-3        # migration rate between neighbours


def build_demography() -> msprime.Demography:
    return msprime.Demography.stepping_stone_model(
        initial_size=[10000] * K,
        migration_rate=M,
    )


def run(random_seed: int = 42) -> dict:
    all_results = {}

    d = build_demography()
    # stepping_stone_model names populations pop_0, pop_1, ...
    pop_scores = {f"pop_{i}": i / (K - 1) for i in range(K)}

    # --- Full gradient: sample from all K populations ---
    print("\n=== Stepping stone: full gradient (all populations sampled) ===")
    samples_full = {f"pop_{i}": N_EACH for i in range(K)}
    out_dir = DATA_DIR / "outputs" / "msprime_stepping_stone_full"
    _, ind_names, pop_labels = simulate(d, samples_full, out_dir,
                                         sequence_length=10e6, random_seed=random_seed)
    pheno = stratified_phenotype(pop_labels, pop_scores, alpha=0.8, random_seed=random_seed)
    lambdas_full = run_experiment(out_dir, ind_names, pop_labels, pheno)
    all_results["full_gradient"] = lambdas_full

    # --- Endpoints only: sample POP_0 and POP_{K-1} only (looks like 2-pop study) ---
    print("\n=== Stepping stone: endpoints only (2-pop illusion) ===")
    n_endpoint = 350
    samples_ends = {"pop_0": n_endpoint, f"pop_{K-1}": n_endpoint}
    out_dir_ends = DATA_DIR / "outputs" / "msprime_stepping_stone_endpoints"
    _, ind_names_e, pop_labels_e = simulate(d, samples_ends, out_dir_ends,
                                             sequence_length=10e6, random_seed=random_seed)
    pheno_e = stratified_phenotype(pop_labels_e, pop_scores, alpha=0.8, random_seed=random_seed)
    lambdas_ends = run_experiment(out_dir_ends, ind_names_e, pop_labels_e, pheno_e)
    all_results["endpoints_only"] = lambdas_ends

    summary_path = DATA_DIR / "outputs" / "msprime_stepping_stone_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n=> Summary: {summary_path}")
    _print_comparison(all_results)
    return all_results


def _print_comparison(results: dict):
    pc_counts = sorted({k for v in results.values() for k in v})
    header = f"{'Scenario':<25}" + "".join(f"  PC{n:>2}" for n in pc_counts)
    print(header)
    for name, lams in results.items():
        row = f"{name:<25}" + "".join(
            f"  {lams.get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        print(row)


if __name__ == "__main__":
    run()
