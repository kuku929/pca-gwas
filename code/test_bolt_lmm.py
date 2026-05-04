"""
BOLT-LMM on all admixture experiments.

Compares GRM-based mixed model correction against PCA correction.
Saves bolt_lambda.json to each experiment's evaluation/ dir.

Usage:
  uv run python test_bolt_lmm.py                  # all experiments
  uv run python test_bolt_lmm.py admixture_gradient three_way_admixture

Install BOLT-LMM once:
  cd /home/krutarth/courses/compbio/project
  wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz
  tar xzf BOLT-LMM_v2.5.tar.gz
"""
import json
import sys
import subprocess
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

OUTPUTS_DIR = Path("/home/krutarth/courses/compbio/data/outputs")
BOLT_BIN    = Path("/home/krutarth/courses/compbio/project/BOLT-LMM_v2.5/bolt")
PREFIX      = "test_chr-2"

ALL_EXPERIMENTS = [
    "admixture_gradient",
    "three_way_admixture",
    "admixture_overwhelm",
    "diff_heritability",
    "admixture_sweep_0pct",
    "admixture_sweep_25pct",
    "admixture_sweep_50pct",
    "admixture_sweep_75pct",
    "admixture_sweep_100pct",
]


def run_experiment(exp: str):
    data_dir = OUTPUTS_DIR / f"test_{exp}"
    eval_dir = data_dir / "evaluation"
    eval_dir.mkdir(exist_ok=True)

    pheno_src = data_dir / "test_chr.pheno1"
    if not pheno_src.exists():
        print(f"[skip] {exp}: test_chr.pheno1 not found")
        return

    # Build FID/IID pheno file
    df = pd.read_csv(pheno_src, sep='\t')
    sample0 = str(df['Sample'].iloc[0])
    if '_' in sample0:
        fid = df['Sample'].str.rsplit('_', n=1).str[0]
        iid = df['Sample'].str.rsplit('_', n=1).str[1]
    else:
        fid = iid = df['Sample']
    liability_col = next(c for c in df.columns if 'liability' in c.lower())
    pheno_file = eval_dir / "bolt_pheno.txt"
    pd.DataFrame({'FID': fid, 'IID': iid, 'PHENO': df[liability_col]}).to_csv(
        pheno_file, sep='\t', index=False
    )

    stats_file = eval_dir / f"bolt_{exp}.stats"
    cmd = [
        str(BOLT_BIN),
        "--bed",               str(data_dir / f"{PREFIX}.bed"),
        "--bim",               str(data_dir / f"{PREFIX}.bim"),
        "--fam",               str(data_dir / f"{PREFIX}.fam"),
        "--phenoFile",         str(pheno_file),
        "--phenoCol",          "PHENO",
        "--lmmForceNonInf",
        "--LDscoresUseChip",
        "--numLeaveOutChunks", "2",
        "--statsFile",         str(stats_file),
    ]
    print(f"\n{'='*55}")
    print(f"=== BOLT-LMM: {exp} ===")
    print(f"{'='*55}")
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[error] BOLT-LMM failed (exit {result.returncode})")
        print(result.stderr[-2000:] if result.stderr else "")
        return

    if not stats_file.exists():
        print(f"[error] stats file not written")
        return

    df_stats = pd.read_csv(stats_file, sep='\t')
    p_col = 'P_BOLT_LMM' if 'P_BOLT_LMM' in df_stats.columns else 'P_BOLT_LMM_INF'
    pvals = pd.to_numeric(df_stats[p_col], errors='coerce').dropna().values
    pvals = pvals[pvals > 0]
    lam = float(np.median(stats.chi2.ppf(1 - pvals, df=1)) / 0.456)

    summary_file = eval_dir / "bolt_lambda.json"
    with open(summary_file, 'w') as f:
        json.dump({"bolt_lmm": round(lam, 4)}, f, indent=2)
    print(f"=> λ = {lam:.4f}  (saved {summary_file.name})")


def main():
    if not BOLT_BIN.exists():
        print(f"[error] BOLT-LMM binary not found at {BOLT_BIN}")
        print("Install:")
        print("  cd /home/krutarth/courses/compbio/project")
        print("  wget https://storage.googleapis.com/broad-alkesgroup-public/BOLT-LMM/downloads/BOLT-LMM_v2.5.tar.gz")
        print("  tar xzf BOLT-LMM_v2.5.tar.gz")
        return

    targets = [a for a in sys.argv[1:] if not a.startswith("--")] or ALL_EXPERIMENTS
    for exp in targets:
        run_experiment(exp)


if __name__ == "__main__":
    main()
