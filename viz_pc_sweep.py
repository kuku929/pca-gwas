#!/usr/bin/env python3
"""
Lambda vs n_PCs line chart across admixture experiments.
Reads pc_sweep_lambdas.json (and bolt_lambda.json if present) from evaluation dirs.

Usage: python viz_pc_sweep.py
"""
import json
from pathlib import Path

import plotly.graph_objects as go

DATA_ROOT = Path("/home/krutarth/courses/compbio/data/outputs")
OUT_FILE  = Path("/home/krutarth/courses/compbio/project/interactive_pc_sweep.html")

EXPERIMENTS = {
    "admixture_gradient":  ("EUR–AFR Admixture Gradient",             "#2166AC"),
    "three_way_admixture": ("Three-Way Admixture",                    "#D6604D"),
    "admixture_overwhelm": ("Admixed Majority Overwhelms",            "#888888"),
    "diff_heritability":   ("Differential Heritability (h²=0.3/0.01)", "#4DAC26"),
}


def main():
    fig = go.Figure()

    # reference line λ = 1
    fig.add_hline(
        y=1.0, line_dash="dash", line_color="black", line_width=1.5,
        annotation_text="λ = 1.0 (no inflation)", annotation_position="top right",
        annotation_font_size=13,
    )

    any_data = False
    for exp, (label, color) in EXPERIMENTS.items():
        json_file = DATA_ROOT / f"test_{exp}" / "evaluation" / "pc_sweep_lambdas.json"
        if not json_file.exists():
            print(f"[skip] {exp}: pc_sweep_lambdas.json not found")
            continue

        data = json.loads(json_file.read_text())
        xs = sorted(int(k) for k in data)
        ys = [data[str(k)] for k in xs]

        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines+markers",
            name=label,
            line=dict(color=color, width=2.5),
            marker=dict(size=10, symbol="circle"),
            hovertemplate=f"<b>{label}</b><br>PCs: %{{x}}<br>λ = %{{y:.4f}}<extra></extra>",
        ))

        # annotate final point with lambda value
        fig.add_annotation(
            x=xs[-1], y=ys[-1],
            text=f"λ={ys[-1]:.3f}",
            showarrow=False,
            xanchor="left", xshift=8,
            font=dict(size=12, color=color),
        )
        any_data = True

        # BOLT-LMM reference point if available (shown as star on the same line)
        bolt_file = DATA_ROOT / f"test_{exp}" / "evaluation" / "bolt_lambda.json"
        if bolt_file.exists():
            bolt_lam = json.loads(bolt_file.read_text()).get("bolt_lmm")
            if bolt_lam is not None:
                fig.add_trace(go.Scatter(
                    x=[xs[-1]], y=[bolt_lam],
                    mode="markers",
                    name=f"{label} — BOLT-LMM",
                    marker=dict(size=16, symbol="star", color=color,
                                line=dict(width=1.5, color="white")),
                    hovertemplate=f"<b>BOLT-LMM ({label})</b><br>λ = {bolt_lam:.4f}<extra></extra>",
                    showlegend=True,
                ))
                fig.add_annotation(
                    x=xs[-1], y=bolt_lam,
                    text=f"BOLT-LMM λ={bolt_lam:.3f}",
                    showarrow=True, arrowhead=2, arrowsize=0.8,
                    ax=30, ay=-30,
                    font=dict(size=11, color=color),
                )

    if not any_data:
        print("No data found. Run test_pc_sweep.py first.")
        return

    fig.update_layout(
        title=(
            "Genomic Inflation Factor (λ) vs Number of PCA Covariates<br>"
            "<sup>PCA correction cannot drive λ→1 for admixed cohorts regardless of PC count. "
            "★ = BOLT-LMM (GRM-based) for comparison.</sup>"
        ),
        xaxis=dict(
            title="Number of PCs used as GWAS covariates",
            tickmode="array",
            tickvals=[0, 2, 5, 10],
            ticktext=["0 PCs<br>(no correction)", "2 PCs", "5 PCs", "10 PCs"],
        ),
        yaxis_title="Genomic inflation factor λ",
        template="plotly_white",
        font=dict(family="Arial", size=14),
        legend=dict(title="Experiment", borderwidth=1, x=0.65, y=0.98),
        width=860, height=580,
        title_font_size=16,
    )

    fig.write_html(OUT_FILE, include_plotlyjs="cdn")
    print(f"wrote {OUT_FILE}")


def plot_admixture_sweep():
    """
    Lambda vs cohort admixture proportion (0→100%).
    Two lines: 0 PCs (no correction) and 10 PCs correction.
    Shows whether/where PCA correction breaks down as cohort composition changes.
    """
    SWEEP_LABELS = ["0pct", "25pct", "50pct", "75pct", "100pct"]
    PCT_VALUES   = [0, 25, 50, 75, 100]

    data = {}
    for label, pct in zip(SWEEP_LABELS, PCT_VALUES):
        json_file = DATA_ROOT / f"test_admixture_sweep_{label}" / "evaluation" / "pc_sweep_lambdas.json"
        if json_file.exists():
            data[pct] = json.loads(json_file.read_text())

    if not data:
        print("[skip] admixture sweep: no pc_sweep_lambdas.json found — run test_admixture_sweep.py first")
        return

    fig = go.Figure()
    fig.add_hline(
        y=1.0, line_dash="dash", line_color="black", line_width=1.5,
        annotation_text="λ = 1.0 (no inflation)", annotation_position="top right",
        annotation_font_size=13,
    )

    PC_LINES = [
        ("0",  "#D6604D", "0 PCs — no correction"),
        ("10", "#2166AC", "10 PCs correction"),
    ]
    for pc_key, color, name in PC_LINES:
        xs = [pct for pct in PCT_VALUES if pct in data and pc_key in data[pct]]
        ys = [data[pct][pc_key] for pct in xs]
        if not xs:
            continue
        fig.add_trace(go.Scatter(
            x=xs, y=ys,
            mode="lines+markers",
            name=name,
            line=dict(color=color, width=2.5),
            marker=dict(size=10),
            hovertemplate="Admixed: %{x}%<br>λ = %{y:.4f}<extra></extra>",
        ))

    fig.update_layout(
        title=(
            "GWAS Inflation (λ) vs Cohort Admixture Proportion<br>"
            "<sup>Fixed N=700. Pure-pop anchor (EUR+AFR) vs increasing admixed fraction. "
            "PCA correction degrades at intermediate admixture; fully admixed cohort "
            "is homogeneous so λ drops back toward 1.</sup>"
        ),
        xaxis=dict(
            title="% admixed individuals in cohort",
            tickmode="array",
            tickvals=PCT_VALUES,
            ticktext=["0%<br>(pure pops)", "25%", "50%", "75%", "100%<br>(fully admixed)"],
        ),
        yaxis_title="Genomic inflation factor λ",
        template="plotly_white",
        font=dict(family="Arial", size=14),
        legend=dict(title="PC correction", borderwidth=1),
        width=860, height=580,
        title_font_size=15,
    )

    out = Path("/home/krutarth/courses/compbio/project/interactive_admixture_sweep.html")
    fig.write_html(out, include_plotlyjs="cdn")
    print(f"wrote {out}")


if __name__ == "__main__":
    import sys
    if "--sweep" in sys.argv:
        plot_admixture_sweep()
    else:
        main()
        plot_admixture_sweep()   # run both by default
