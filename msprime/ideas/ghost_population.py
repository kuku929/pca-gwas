"""
Idea 2: Ghost population — unsampled ancestry source.

Ghost:   A + ADMIX_AC. C never sampled → PCA only sees A variation.
         C-ancestry in ADMIX_AC is orthogonal to all PC axes → λ uncorrectable.
Control: A + C + ADMIX_AC. C sampled → PCA can see A–C axis → λ → 1.
"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent.parent))

import json
import msprime
from simulate import simulate
from phenotype import stratified_phenotype
from gwas import run_experiment

PROJECT_ROOT = Path(__file__).parent.parent.parent
DATA_DIR = PROJECT_ROOT / "data"


def build_demography():
    d = msprime.Demography()
    d.add_population(name="A",        initial_size=10000)
    d.add_population(name="C",        initial_size=10000)
    d.add_population(name="ADMIX_AC", initial_size=10000)
    d.add_population(name="ANC",      initial_size=20000)
    d.add_admixture(time=50, derived="ADMIX_AC", ancestral=["A", "C"],
                    proportions=[0.5, 0.5])
    d.add_population_split(time=1000, derived=["A", "C"], ancestral="ANC")
    return d


def run(random_seed: int = 42) -> dict:
    all_results = {}
    d = build_demography()

    # --- Ghost: A + ADMIX_AC only (C unsampled) ---
    print("\n=== Ghost (A + ADMIX_AC, C unsampled) ===")
    out_dir = DATA_DIR / "outputs" / "msprime_ghost"
    _, ind_names, pop_labels = simulate(d, {"A": 400, "ADMIX_AC": 300}, out_dir,
                                        sequence_length=10e6, random_seed=random_seed)
    pop_scores = {"A": 0.0, "ADMIX_AC": 1.0}
    pheno = stratified_phenotype(pop_labels, pop_scores, alpha=0.8, random_seed=random_seed)
    all_results["ghost"] = run_experiment(out_dir, ind_names, pop_labels, pheno)

    # --- Control: A + C + ADMIX_AC (C sampled) ---
    print("\n=== Control (A + C + ADMIX_AC, C sampled) ===")
    out_dir_ctrl = DATA_DIR / "outputs" / "msprime_ghost_control"
    _, ind_names_c, pop_labels_c = simulate(d, {"A": 400, "C": 150, "ADMIX_AC": 300},
                                             out_dir_ctrl, sequence_length=10e6,
                                             random_seed=random_seed)
    pop_scores_ctrl = {"A": 0.0, "C": 1.0, "ADMIX_AC": 1.0}
    pheno_c = stratified_phenotype(pop_labels_c, pop_scores_ctrl, alpha=0.8,
                                   random_seed=random_seed)
    all_results["control"] = run_experiment(out_dir_ctrl, ind_names_c, pop_labels_c, pheno_c)

    summary_path = DATA_DIR / "outputs" / "msprime_ghost_summary.json"
    with open(summary_path, "w") as f:
        json.dump(all_results, f, indent=2)
    print(f"\n=> Summary: {summary_path}")
    _print_comparison(all_results)
    return all_results


def _print_comparison(results: dict):
    pc_counts = sorted({k for v in results.values() for k in v if k > 0})
    header = f"{'Scenario':<20}" + "".join(f"  PC{n:>2}" for n in pc_counts)
    print(header)
    for name, lams in results.items():
        row = f"{name:<20}" + "".join(
            f"  {lams.get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        print(row)


if __name__ == "__main__":
    run()
