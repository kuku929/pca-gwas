"""Core simulation: demographic model → VCF + individual metadata."""
import msprime
import numpy as np
from pathlib import Path


def simulate(
    demography: msprime.Demography,
    samples: dict,
    output_path: Path,
    sequence_length: float = 10e6,
    recombination_rate: float = 1e-8,
    mutation_rate: float = 1.29e-8,
    random_seed: int = 42,
) -> tuple:
    """
    Run sim_ancestry + sim_mutations. Write VCF to output_path/sim.vcf.

    Returns (mts, individual_names, pop_labels):
      individual_names: ["ind_0", "ind_1", ...] matching VCF columns
      pop_labels:       population name for each individual
    """
    output_path = Path(output_path)
    output_path.mkdir(parents=True, exist_ok=True)

    ts = msprime.sim_ancestry(
        samples=samples,
        demography=demography,
        sequence_length=sequence_length,
        recombination_rate=recombination_rate,
        random_seed=random_seed,
    )
    mts = msprime.sim_mutations(ts, rate=mutation_rate, random_seed=random_seed + 1)

    pop_id_to_name = {pop.id: pop.metadata["name"] for pop in mts.populations()}
    individual_names = []
    pop_labels = []
    for i, ind in enumerate(mts.individuals()):
        individual_names.append(f"ind_{i}")
        pop_id = mts.node(ind.nodes[0]).population
        pop_labels.append(pop_id_to_name[pop_id])

    vcf_path = output_path / "sim.vcf"
    with open(vcf_path, "w") as f:
        mts.write_vcf(f, individual_names=individual_names)

    n_snps = mts.num_sites
    print(f"Simulated {mts.num_individuals} individuals, {n_snps} sites → {vcf_path}")
    return mts, individual_names, pop_labels
