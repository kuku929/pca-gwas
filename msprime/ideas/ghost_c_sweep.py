"""
Ghost C-sample sweep.

Fix A(400) + ADMIX_AC(300). Sweep number of pure C samples added: 0, 25, 50, 100, 200, 400.
Phenotype: A=0, C=1, ADMIX_AC=1.
Shows whether adding pure C samples helps PCA correct stratification.
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

C_COUNTS = [0, 25, 50, 100, 200, 400]
N_A = 400
N_ADMIX = 300


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
    pop_scores = {"A": 0.0, "C": 1.0, "ADMIX_AC": 1.0}

    for n_c in C_COUNTS:
        print(f"\n=== C samples = {n_c} (total = {N_A + N_ADMIX + n_c}) ===")
        out_dir = DATA_DIR / "outputs" / f"msprime_ghost_c{n_c}"
        samples = {"A": N_A, "ADMIX_AC": N_ADMIX}
        if n_c > 0:
            samples["C"] = n_c
        _, ind_names, pop_labels = simulate(d, samples, out_dir,
                                            sequence_length=10e6, random_seed=random_seed)
        pheno = stratified_phenotype(pop_labels, pop_scores, alpha=0.8,
                                     random_seed=random_seed)
        lambdas = run_experiment(out_dir, ind_names, pop_labels, pheno, run_bolt_lmm=True)
        all_results[n_c] = lambdas

    summary_path = DATA_DIR / "outputs" / "msprime_ghost_c_sweep_summary.json"
    with open(summary_path, "w") as f:
        json.dump({str(k): v for k, v in all_results.items()}, f, indent=2)
    print(f"\n=> Summary: {summary_path}")
    _print_table(all_results)
    return all_results


def _print_table(results: dict):
    pc_counts = sorted({k for v in results.values() for k in v if isinstance(k, int) and k > 0})
    print(f"\n{'C samples':>10}" + "".join(f"  PC{n:>2}" for n in pc_counts))
    for n_c in sorted(results):
        lams = results[n_c]
        row = f"{n_c:>10}" + "".join(
            f"  {lams.get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        print(row)


if __name__ == "__main__":
    run()
