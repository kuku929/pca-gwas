#!/usr/bin/env python3
"""
Interactive PCA visualizations for HAPNEST admixture experiments.
Writes Plotly HTML files to each experiment's evaluation/ directory.

Usage:
    python viz_pca.py                          # all three experiments
    python viz_pca.py admixture_gradient       # one experiment
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import plotly.express as px
import plotly.graph_objects as go
from plotly.subplots import make_subplots

DATA_ROOT = Path("/home/krutarth/courses/compbio/data/outputs")
CHR = "2"

EXPERIMENTS = ["admixture_gradient", "three_way_admixture", "admixture_overwhelm"]

# ── Per-experiment color palettes ─────────────────────────────────────────────
# Gradient: RdBu diverging — encodes EUR→AFR spectrum visually
COLORS = {
    "admixture_gradient": {
        "pure_EUR":    "#2166AC",
        "admix_75EUR": "#74ADD1",
        "admix_50EUR": "#9B6BB5",
        "admix_25EUR": "#F4A582",
        "pure_AFR":    "#D6604D",
    },
    "three_way_admixture": {
        "pure_EUR":       "#2166AC",
        "pure_AFR":       "#D6604D",
        "pure_EAS":       "#1A9850",
        "admix_EUR_AFR":  "#762A83",
        "admix_EUR_EAS":  "#4393C3",
        "admix_AFR_EAS":  "#A6D96A",
        "admix_three_way":"#FEE090",
    },
    "admixture_overwhelm": {
        "pure_EUR":         "#2166AC",
        "pure_AFR":         "#D6604D",
        "admixed_majority": "#BBBBBB",
    },
}

# Gradient only: known AFR ancestry % per group (ground truth)
AFR_PROPORTION = {
    "pure_EUR":    0,
    "admix_75EUR": 25,
    "admix_50EUR": 50,
    "admix_25EUR": 75,
    "pure_AFR":    100,
}

TITLES = {
    "admixture_gradient": "EUR–AFR Admixture Gradient",
    "three_way_admixture": "Three-Way Admixture (EUR · AFR · EAS)",
    "admixture_overwhelm": "Admixed Majority Overwhelms PCA",
}


# ── Data loading ──────────────────────────────────────────────────────────────

def load_experiment(exp_name: str):
    exp_dir  = DATA_ROOT / f"test_{exp_name}"
    eval_dir = exp_dir / "evaluation"
    prefix   = f"test_chr-{CHR}"

    pca = pd.read_csv(eval_dir / f"{prefix}.syn.pca.eigenvec", sep=r"\s+")
    # plink2 names first col '#FID'
    pca.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, len(pca.columns) - 1)]

    eigenval = pd.read_csv(
        eval_dir / f"{prefix}.syn.pca.eigenval", header=None, names=["val"]
    )["val"]

    labels = pd.read_csv(
        exp_dir / f"{prefix}.sample", header=None, names=["Population"]
    )
    pca["Population"] = labels["Population"].values
    return pca, eigenval


# ── Plot: PC1 vs PC2 scatter ──────────────────────────────────────────────────

def plot_pca_scatter(pca: pd.DataFrame, exp_name: str, eval_dir: Path):
    colors = COLORS.get(exp_name, {})
    title  = TITLES.get(exp_name, exp_name)

    fig = px.scatter(
        pca, x="PC1", y="PC2",
        color="Population",
        color_discrete_map=colors,
        hover_data={"IID": True, "PC1": ":.4f", "PC2": ":.4f", "Population": True},
        title=f"PCA — {title}",
        template="plotly_white",
        labels={"PC1": "PC 1", "PC2": "PC 2"},
    )
    fig.update_traces(
        marker=dict(size=7, opacity=0.85, line=dict(width=0.5, color="white"))
    )
    fig.update_layout(
        font=dict(family="Arial", size=14),
        legend=dict(title="Population group", itemsizing="constant", borderwidth=1),
        width=820, height=660,
        title_font_size=17,
    )
    out = eval_dir / f"interactive_pca_scatter_{exp_name}.html"
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  {out.name}")


# ── Plot: scree (variance explained) ─────────────────────────────────────────

def plot_scree(eigenval: pd.Series, exp_name: str, eval_dir: Path):
    pct  = eigenval / eigenval.sum() * 100
    cum  = pct.cumsum()
    labs = [f"PC{i+1}" for i in range(len(pct))]

    fig = make_subplots(specs=[[{"secondary_y": True}]])
    fig.add_trace(
        go.Bar(x=labs, y=pct.values, name="% variance",
               marker_color="#2166AC", opacity=0.75),
        secondary_y=False,
    )
    fig.add_trace(
        go.Scatter(x=labs, y=cum.values, name="cumulative %",
                   mode="lines+markers",
                   line=dict(color="#D6604D", width=2),
                   marker=dict(size=6)),
        secondary_y=True,
    )
    fig.update_layout(
        title=f"Variance Explained — {TITLES.get(exp_name, exp_name)}",
        template="plotly_white",
        font=dict(family="Arial", size=14),
        legend=dict(x=0.65, y=0.45, borderwidth=1),
        width=760, height=500,
        title_font_size=17,
    )
    fig.update_yaxes(title_text="% variance per PC", secondary_y=False)
    fig.update_yaxes(title_text="Cumulative %", secondary_y=True, range=[0, 108])
    out = eval_dir / f"interactive_scree_{exp_name}.html"
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  {out.name}")


# ── Plot: PC1 as ancestry proxy (gradient only) ───────────────────────────────

def plot_ancestry_proxy(pca: pd.DataFrame, eval_dir: Path):
    df = pca.copy()
    df["AFR_pct"] = df["Population"].map(AFR_PROPORTION)
    # Ordered x-axis labels
    label_map = {k: f"{v}% AFR  ({k})" for k, v in AFR_PROPORTION.items()}
    df["Group"] = df["Population"].map(label_map)
    order = [label_map[k] for k in sorted(AFR_PROPORTION, key=AFR_PROPORTION.get)]

    fig = px.violin(
        df, x="Group", y="PC1",
        color="Population",
        color_discrete_map=COLORS["admixture_gradient"],
        category_orders={"Group": order},
        box=True, points="all",
        hover_data=["IID"],
        title=(
            "PC 1 Score vs Known Admixture Proportion<br>"
            "<sup>A linear relationship proves PC1 measures ancestry, not independent biological signal</sup>"
        ),
        template="plotly_white",
        labels={"PC1": "PC 1 score", "Group": "Known AFR ancestry proportion"},
    )
    fig.update_traces(marker=dict(size=3, opacity=0.45), width=0.75)
    fig.update_layout(
        font=dict(family="Arial", size=13),
        showlegend=False,
        width=900, height=600,
        title_font_size=16,
    )
    out = eval_dir / "interactive_ancestry_proxy.html"
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  {out.name}")


# ── Plot: 3D scatter (three-way only) ─────────────────────────────────────────

def plot_pca_3d(pca: pd.DataFrame, eval_dir: Path):
    colors = COLORS["three_way_admixture"]

    fig = px.scatter_3d(
        pca, x="PC1", y="PC2", z="PC3",
        color="Population",
        color_discrete_map=colors,
        hover_data={"IID": True},
        title=(
            "3-Way Admixture — PCA in 3D<br>"
            "<sup>Admixed groups fill triangle edges and interior — PCA axes capture ancestry geometry</sup>"
        ),
        template="plotly_white",
        opacity=0.82,
    )
    fig.update_traces(marker=dict(size=4))
    fig.update_layout(
        font=dict(family="Arial", size=13),
        legend=dict(title="Population group", itemsizing="constant", borderwidth=1),
        width=880, height=720,
        title_font_size=16,
        scene=dict(
            xaxis_title="PC 1",
            yaxis_title="PC 2",
            zaxis_title="PC 3",
            bgcolor="white",
        ),
    )
    out = eval_dir / "interactive_pca_3d_three_way.html"
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"  {out.name}")


# ── Runner ────────────────────────────────────────────────────────────────────

def run_experiment(exp_name: str):
    exp_dir  = DATA_ROOT / f"test_{exp_name}"
    eval_dir = exp_dir / "evaluation"
    sample   = exp_dir / f"test_chr-{CHR}.sample"

    if not sample.exists():
        print(f"[skip] {exp_name}: {sample} not found — run the test first")
        return

    print(f"\n=== {exp_name} ===")
    pca, eigenval = load_experiment(exp_name)

    plot_pca_scatter(pca, exp_name, eval_dir)
    plot_scree(eigenval, exp_name, eval_dir)

    if exp_name == "admixture_gradient":
        plot_ancestry_proxy(pca, eval_dir)

    if exp_name == "three_way_admixture":
        plot_pca_3d(pca, eval_dir)


if __name__ == "__main__":
    targets = sys.argv[1:] if len(sys.argv) > 1 else EXPERIMENTS
    for exp in targets:
        run_experiment(exp)
