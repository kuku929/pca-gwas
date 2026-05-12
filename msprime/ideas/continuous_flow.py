"""
Idea 4: Pulse admixture vs continuous gene flow.

Both scenarios have EUR + AFR + ADMIX cohort with the same *mean* ancestry proportion.
The difference: where does the variance in ancestry come from?

Scenario A (pulse): Single event at T=50 gens. ADMIX individuals have exactly 50% EUR / 50% AFR.
  PCA sees 3 tight clusters. 2 PCs largely correct.

Scenario B (continuous): EUR↔AFR symmetric migration m=0.02/gen for T=50 gens, then stopped.
  EUR and AFR individuals have *variable* amounts of cross-ancestry (distribution, not fixed).
  PCA sees 2 fuzzy clusters. λ remains elevated because variance in ancestry ≠ mean ancestry.
  PCA axis captures mean ancestry difference; variance around that axis is uncorrected.
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

MIGRATION_RATE = 0.02   # per generation, symmetric EUR↔AFR
MIGRATION_GENS = 50     # gens of continuous flow (same as pulse time T)


def build_pulse_demography():
    d = msprime.Demography()
    d.add_population(name="EUR",   initial_size=10000)
    d.add_population(name="AFR",   initial_size=14000)
    d.add_population(name="ADMIX", initial_size=10000)
    d.add_population(name="ANC",   initial_size=20000)
    d.add_admixture(time=MIGRATION_GENS, derived="ADMIX",
                    ancestral=["EUR", "AFR"], proportions=[0.5, 0.5])
    d.add_population_split(time=2000, derived=["EUR", "AFR"], ancestral="ANC")
    return d


def build_continuous_demography():
    d = msprime.Demography()
    d.add_population(name="EUR", initial_size=10000)
    d.add_population(name="AFR", initial_size=14000)
    d.add_population(name="ANC", initial_size=20000)
    # Migration active from present back to T=MIGRATION_GENS (backwards-in-time)
    d.set_symmetric_migration_rate(["EUR", "AFR"], MIGRATION_RATE)
    # Turn off migration before MIGRATION_GENS (migration "started" MIGRATION_GENS ago)
    d.add_symmetric_migration_rate_change(
        time=MIGRATION_GENS, populations=["EUR", "AFR"], rate=0.0
    )
    d.add_population_split(time=2000, derived=["EUR", "AFR"], ancestral="ANC")
    return d


def run(random_seed: int = 42) -> dict:
    all_results = {}

    # --- Pulse scenario ---
    print(f"\n=== Pulse admixture (T={MIGRATION_GENS} gens, 50/50 EUR+AFR) ===")
    out_dir_pulse = DATA_DIR / "outputs" / "msprime_flow_pulse"
    d_pulse = build_pulse_demography()
    samples_pulse = {"EUR": 200, "AFR": 200, "ADMIX": 300}
    _, ind_names_p, pop_labels_p = simulate(d_pulse, samples_pulse, out_dir_pulse,
                                             sequence_length=10e6, random_seed=random_seed)
    pop_scores_pulse = {"EUR": 1.0, "AFR": 0.0, "ADMIX": 0.5}
    pheno_p = stratified_phenotype(pop_labels_p, pop_scores_pulse, alpha=0.8,
                                    random_seed=random_seed)
    lambdas_pulse = run_experiment(out_dir_pulse, ind_names_p, pop_labels_p, pheno_p)
    all_results["pulse"] = lambdas_pulse

    # --- Continuous flow scenario ---
    print(f"\n=== Continuous gene flow (m={MIGRATION_RATE}/gen for {MIGRATION_GENS} gens) ===")
    out_dir_cont = DATA_DIR / "outputs" / "msprime_flow_continuous"
    d_cont = build_continuous_demography()
    # No separate ADMIX population — EUR and AFR labels are fuzzy after gene flow
    samples_cont = {"EUR": 350, "AFR": 350}
    _, ind_names_c, pop_labels_c = simulate(d_cont, samples_cont, out_dir_cont,
                                             sequence_length=10e6, random_seed=random_seed)
    pop_scores_cont = {"EUR": 1.0, "AFR": 0.0}
    pheno_c = stratified_phenotype(pop_labels_c, pop_scores_cont, alpha=0.8,
                                    random_seed=random_seed)
    lambdas_cont = run_experiment(out_dir_cont, ind_names_c, pop_labels_c, pheno_c)
    all_results["continuous"] = lambdas_cont

    summary_path = DATA_DIR / "outputs" / "msprime_flow_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n=> Summary: {summary_path}")
    _print_comparison(all_results)
    return all_results


def _print_comparison(results: dict):
    pc_counts = sorted({k for v in results.values() for k in v})
    header = f"{'Scenario':<15}" + "".join(f"  PC{n:>2}" for n in pc_counts)
    print(header)
    for name, lams in results.items():
        row = f"{name:<15}" + "".join(
            f"  {lams.get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        print(row)


if __name__ == "__main__":
    run()
