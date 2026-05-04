#!/usr/bin/env python3
"""
Interactive QQ plots for GWAS results from HAPNEST admixture experiments.
Handles both plink2 (.glm.linear) and BOLT-LMM (.stats) output formats.

Usage:
    python viz_gwas.py                          # all three experiments (plink2)
    python viz_gwas.py admixture_gradient       # one experiment
    python viz_gwas.py --bolt admixture_gradient  # BOLT-LMM QQ for gradient
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.graph_objects as go
from scipy import stats

DATA_ROOT = Path("/home/krutarth/courses/compbio/data/outputs")
CHR = "2"

EXPERIMENTS = ["admixture_gradient", "three_way_admixture", "admixture_overwhelm"]

TITLES = {
    "admixture_gradient": "EUR–AFR Admixture Gradient",
    "three_way_admixture": "Three-Way Admixture (EUR · AFR · EAS)",
    "admixture_overwhelm": "Admixed Majority Overwhelms PCA",
}

# plink2 GWAS output extensions to search for
GWAS_EXTENSIONS = [
    ".glm.linear",
    ".glm.logistic",
    ".glm.logistic.hybrid",
    ".assoc.linear",
    ".assoc.logistic",
]


def find_gwas_file(exp_dir: Path, prefix: str) -> Path | None:
    eval_dir = exp_dir / "evaluation"

    # Prefer pc_sweep output at highest PC count (most corrected, best comparison)
    for n in [10, 5, 2, 0]:
        for pat in [f"{prefix}.gwas_pc{n}.PHENO.glm.linear",
                    f"{prefix}.gwas_pc{n}.PHENO1.glm.linear"]:
            f = eval_dir / pat
            if f.exists():
                return f

    # Fall back to validate's plink2 output (trait1 naming)
    for ext in GWAS_EXTENSIONS:
        for search_dir in [eval_dir, exp_dir]:
            candidate = search_dir / f"{prefix}{ext}"
            if candidate.exists():
                return candidate

    # Last resort: any .glm.* file in eval dir
    hits = sorted(eval_dir.glob("*.glm.*"))
    return hits[0] if hits else None


def compute_lambda(p_values: np.ndarray) -> float:
    """Genomic inflation factor λ = median(χ²_obs) / 0.456"""
    chi2_obs = stats.chi2.ppf(1.0 - p_values, df=1)
    return float(np.median(chi2_obs) / 0.456)


def plot_qq(p_values: np.ndarray, exp_name: str, eval_dir: Path):
    p_values = p_values[~np.isnan(p_values) & (p_values > 0)]
    n   = len(p_values)
    lam = compute_lambda(p_values)

    obs = -np.log10(np.sort(p_values))                  # ascending p → descending -log10
    exp = -np.log10(np.arange(1, n + 1) / (n + 1))     # expected uniform quantiles

    maxval = max(float(obs.max()), float(exp.max())) * 1.08

    fig = go.Figure()

    # Null diagonal
    fig.add_trace(go.Scatter(
        x=[0, maxval], y=[0, maxval],
        mode="lines",
        line=dict(color="#D6604D", dash="dash", width=1.5),
        name="Expected (null)",
        hoverinfo="skip",
    ))

    # 95% confidence band (Kolmogorov–Smirnov envelope)
    ranks  = np.arange(1, n + 1)
    ci_lo  = -np.log10(stats.beta.ppf(0.975, ranks, n - ranks + 1))
    ci_hi  = -np.log10(stats.beta.ppf(0.025, ranks, n - ranks + 1))
    x_band = exp[::-1]  # exp is descending; reverse to match ci arrays

    fig.add_trace(go.Scatter(
        x=np.concatenate([x_band, x_band[::-1]]),
        y=np.concatenate([ci_lo,  ci_hi[::-1]]),
        fill="toself",
        fillcolor="rgba(150,150,150,0.20)",
        line=dict(width=0),
        name="95% CI",
        hoverinfo="skip",
    ))

    # Observed points
    fig.add_trace(go.Scatter(
        x=exp, y=obs,
        mode="markers",
        marker=dict(size=4, color="#2166AC", opacity=0.7),
        name="Observed",
        text=[f"p = {p:.2e}" for p in np.sort(p_values)],
        hovertemplate="Expected −log₁₀(p): %{x:.3f}<br>Observed −log₁₀(p): %{y:.3f}<br>%{text}",
    ))

    # Lambda annotation box
    lam_color = "#D6604D" if lam > 1.05 else "#1A9850"
    fig.add_annotation(
        x=0.03, y=0.97, xref="paper", yref="paper",
        text=f"λ = {lam:.3f}",
        showarrow=False,
        font=dict(size=16, color=lam_color, family="Arial Black"),
        bgcolor="white",
        bordercolor=lam_color,
        borderwidth=1.5,
        borderpad=6,
    )

    fig.update_layout(
        title=(
            f"QQ Plot — {TITLES.get(exp_name, exp_name)}<br>"
            f"<sup>λ &gt; 1.0 = GWAS inflation = PCA failed to correct population stratification"
            f"  (ground truth heritability = 0.1 uniform across all groups)</sup>"
        ),
        xaxis_title="Expected −log₁₀(p)",
        yaxis_title="Observed −log₁₀(p)",
        template="plotly_white",
        font=dict(family="Arial", size=14),
        legend=dict(x=0.02, y=0.72, borderwidth=1),
        width=720, height=660,
        title_font_size=16,
    )

    out = eval_dir / f"interactive_qqplot_{exp_name}.html"
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  {out.name}  (λ = {lam:.4f})")


def run_experiment(exp_name: str):
    exp_dir  = DATA_ROOT / f"test_{exp_name}"
    eval_dir = exp_dir / "evaluation"
    prefix   = f"test_chr-{CHR}"

    gwas_file = find_gwas_file(exp_dir, prefix)
    if gwas_file is None:
        print(f"[skip] {exp_name}: no GWAS file found — run generate_pheno() + validate() first")
        return

    print(f"\n=== {exp_name}  ({gwas_file.name}) ===")

    df = pd.read_csv(gwas_file, sep="\t")
    df.columns = [c.lstrip("#") for c in df.columns]

    if "TEST" in df.columns:
        df = df[df["TEST"] == "ADD"]

    p_col = "P" if "P" in df.columns else df.columns[-1]
    p_values = pd.to_numeric(df[p_col], errors="coerce").dropna().values

    plot_qq(p_values, exp_name, eval_dir)


def run_bolt_experiment(exp_name: str):
    """QQ plot from BOLT-LMM .stats output."""
    eval_dir   = DATA_ROOT / f"test_{exp_name}" / "evaluation"
    stats_file = eval_dir / f"bolt_{exp_name}.stats"

    if not stats_file.exists():
        print(f"[skip] {exp_name}: bolt_{exp_name}.stats not found — run test_bolt_lmm.py first")
        return

    print(f"\n=== {exp_name} [BOLT-LMM] ===")
    df = pd.read_csv(stats_file, sep="\t")
    p_col = "P_BOLT_LMM" if "P_BOLT_LMM" in df.columns else "P_BOLT_LMM_INF"
    p_values = pd.to_numeric(df[p_col], errors="coerce").dropna().values
    plot_qq(p_values, exp_name + "_bolt", eval_dir)


if __name__ == "__main__":
    args = sys.argv[1:]
    bolt_mode = "--bolt" in args
    targets = [a for a in args if not a.startswith("--")] or EXPERIMENTS
    for exp in targets:
        if bolt_mode:
            run_bolt_experiment(exp)
        else:
            run_experiment(exp)
