"""Phenotype assignment: ancestry-correlated trait (pure stratification, no true genetic effects)."""
import numpy as np
import pandas as pd
from pathlib import Path


def stratified_phenotype(
    pop_labels: list[str],
    pop_scores: dict[str, float],
    alpha: float = 0.8,
    random_seed: int = 0,
) -> np.ndarray:
    """
    Phenotype = alpha * ancestry_score + sqrt(1-alpha^2) * noise.
    alpha controls confounding strength (0=none, 1=fully determined by ancestry).
    pop_scores maps population name → numeric ancestry score.
    """
    rng = np.random.default_rng(random_seed)
    scores = np.array([pop_scores[p] for p in pop_labels], dtype=float)
    # Standardize scores to zero mean, unit variance
    scores = (scores - scores.mean()) / (scores.std() + 1e-9)
    noise = rng.standard_normal(len(pop_labels))
    return alpha * scores + np.sqrt(max(0, 1 - alpha ** 2)) * noise


def write_pheno(output_dir: Path, individual_names: list[str], phenotype: np.ndarray) -> Path:
    """Write plink2-format phenotype file (#FID IID PHENO)."""
    pheno_path = output_dir / "pheno.txt"
    pd.DataFrame({
        "#FID": individual_names,
        "IID": individual_names,
        "PHENO": phenotype,
    }).to_csv(pheno_path, sep="\t", index=False)
    return pheno_path
