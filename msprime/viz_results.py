"""
Visualize msprime PCA failure experiment results.
Reads lambdas.json files from data/outputs/msprime_*/
Produces: HTML plots + printed tables.
"""
import json
import sys
from pathlib import Path
import plotly.graph_objects as go
from plotly.subplots import make_subplots

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data" / "outputs"
OUT_DIR = DATA_DIR / "msprime_plots"
OUT_DIR.mkdir(exist_ok=True)

PC_COUNTS = [2, 5, 10, 20]
COLORS = ["#e41a1c", "#377eb8", "#4daf4a", "#984ea3", "#ff7f00", "#a65628", "#f781bf"]


def load_lambdas(exp_dir: Path) -> dict:
    p = exp_dir / "lambdas.json"
    if not p.exists():
        return {}
    with open(p) as f:
        return {int(k): v for k, v in json.load(f).items()}


# ── Experiment 1: Admixture Timing ──────────────────────────────────────────

def plot_timing():
    summary = DATA_DIR.parent / "outputs" / "msprime_timing_summary.json"
    if summary.exists():
        with open(summary) as f:
            raw = json.load(f)
        data = {int(T): {int(k): v for k, v in lams.items()} for T, lams in raw.items()}
    else:
        # Fallback: read per-experiment dirs
        data = {}
        for d in sorted(DATA_DIR.glob("msprime_timing_T*")):
            T = int(d.name.split("T")[-1])
            lams = load_lambdas(d)
            if lams:
                data[T] = lams

    if not data:
        print("  timing: no results yet")
        return

    T_vals = sorted(data)
    fig = go.Figure()
    for i, n in enumerate(PC_COUNTS):
        y = [data[T].get(n) for T in T_vals]
        if any(v is not None for v in y):
            fig.add_trace(go.Scatter(
                x=T_vals, y=y, mode="lines+markers",
                name=f"PC{n}", line=dict(color=COLORS[i]),
                marker=dict(size=8),
            ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="λ=1 (null)")
    fig.update_layout(
        title="Admixture Timing Sweep: λ vs Generations Since Admixture",
        xaxis=dict(
            title="Admixture time T (generations ago)",
            type="log",
            tickmode="array",
            tickvals=T_vals,
            ticktext=[str(T) for T in T_vals],
        ),
        yaxis_title="Genomic inflation factor λ",
        legend_title="PC covariates",
        template="plotly_white",
        width=800, height=500,
    )
    out = OUT_DIR / "timing_lambda.html"
    fig.write_html(out)
    print(f"  [timing] saved {out}")
    _print_table("Admixture Timing", {str(T): data[T] for T in T_vals}, row_label="T (gens)")


# ── Experiment 2: Ghost Population ──────────────────────────────────────────

def plot_ghost():
    summary = DATA_DIR.parent / "outputs" / "msprime_ghost_summary.json"
    if summary.exists():
        with open(summary) as f:
            data = {k: {int(n): v for n, v in lams.items()} for k, lams in json.load(f).items()}
    else:
        data = {}
        for key, dirname in [("ghost", "msprime_ghost"), ("ghost_control", "msprime_ghost_control")]:
            lams = load_lambdas(DATA_DIR / dirname)
            if lams:
                data[key] = lams

    if not data:
        print("  ghost: no results yet")
        return

    fig = go.Figure()
    styles = {"ghost": dict(color="#e41a1c", dash="solid"), "ghost_control": dict(color="#377eb8", dash="dash")}
    labels = {"ghost": "Ghost (C unsampled)", "ghost_control": "Control (C sampled)"}
    for key, lams in data.items():
        pcs = [n for n in sorted(lams) if n > 0]
        fig.add_trace(go.Scatter(
            x=pcs, y=[lams[n] for n in pcs], mode="lines+markers",
            name=labels.get(key, key), line=styles.get(key, {}),
            marker=dict(size=8),
        ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="λ=1 (null)")
    fig.update_layout(
        title="Ghost Population: λ vs PC Count<br><sup>Ghost=C unsampled (PCA blind spot) vs Control=C sampled</sup>",
        xaxis_title="Number of PC covariates",
        yaxis_title="Genomic inflation factor λ",
        legend_title="Scenario",
        template="plotly_white",
        width=800, height=500,
    )
    out = OUT_DIR / "ghost_lambda.html"
    fig.write_html(out)
    print(f"  [ghost] saved {out}")
    _print_table("Ghost Population", data)


# ── Experiment 3: Stepping Stone ────────────────────────────────────────────

def plot_stepping():
    summary = DATA_DIR.parent / "outputs" / "msprime_stepping_stone_summary.json"
    if summary.exists():
        with open(summary) as f:
            data = {k: {int(n): v for n, v in lams.items()} for k, lams in json.load(f).items()}
    else:
        data = {}
        for key, dirname in [("full_gradient", "msprime_stepping_stone_full"),
                              ("endpoints_only", "msprime_stepping_stone_endpoints")]:
            lams = load_lambdas(DATA_DIR / dirname)
            if lams:
                data[key] = lams

    if not data:
        print("  stepping_stone: no results yet")
        return

    fig = go.Figure()
    styles = {
        "full_gradient":   dict(color="#e41a1c", dash="solid"),
        "endpoints_only":  dict(color="#377eb8", dash="dash"),
    }
    labels = {
        "full_gradient":  "Full gradient (all 16 pops sampled)",
        "endpoints_only": "Endpoints only (2-pop illusion)",
    }
    for key, lams in data.items():
        pcs = [n for n in sorted(lams) if n > 0]
        fig.add_trace(go.Scatter(
            x=pcs, y=[lams[n] for n in pcs], mode="lines+markers",
            name=labels.get(key, key), line=styles.get(key, {}),
            marker=dict(size=8),
        ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="λ=1 (null)")
    fig.update_layout(
        title="Isolation by Distance (Stepping Stone, K=16):<br>λ vs PC Count",
        xaxis_title="Number of PC covariates",
        yaxis_title="Genomic inflation factor λ",
        legend_title="Sampling strategy",
        template="plotly_white",
        width=800, height=500,
    )
    out = OUT_DIR / "stepping_stone_lambda.html"
    fig.write_html(out)
    print(f"  [stepping] saved {out}")
    _print_table("Stepping Stone", data)


# ── Experiment 4: Continuous Flow ───────────────────────────────────────────

def plot_flow():
    summary = DATA_DIR.parent / "outputs" / "msprime_flow_summary.json"
    if summary.exists():
        with open(summary) as f:
            data = {k: {int(n): v for n, v in lams.items()} for k, lams in json.load(f).items()}
    else:
        data = {}
        for key, dirname in [("pulse", "msprime_flow_pulse"), ("continuous", "msprime_flow_continuous")]:
            lams = load_lambdas(DATA_DIR / dirname)
            if lams:
                data[key] = lams

    if not data:
        print("  flow: no results yet")
        return

    fig = go.Figure()
    styles = {
        "pulse":      dict(color="#4daf4a", dash="solid"),
        "continuous": dict(color="#984ea3", dash="dash"),
    }
    labels = {
        "pulse":      "Pulse admixture (T=50, 50/50 EUR+AFR)",
        "continuous": f"Continuous gene flow (m=0.02/gen, 50 gens)",
    }
    for key, lams in data.items():
        pcs = [n for n in sorted(lams) if n > 0]
        fig.add_trace(go.Scatter(
            x=pcs, y=[lams[n] for n in pcs], mode="lines+markers",
            name=labels.get(key, key), line=styles.get(key, {}),
            marker=dict(size=8),
        ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="λ=1 (null)")
    fig.update_layout(
        title="Pulse vs Continuous Gene Flow: λ vs PC Count",
        xaxis_title="Number of PC covariates",
        yaxis_title="Genomic inflation factor λ",
        legend_title="Admixture mechanism",
        template="plotly_white",
        width=800, height=500,
    )
    out = OUT_DIR / "flow_lambda.html"
    fig.write_html(out)
    print(f"  [flow] saved {out}")
    _print_table("Pulse vs Continuous Flow", data)


# ── Combined Summary ─────────────────────────────────────────────────────────

def plot_summary():
    """Single figure showing PC10 λ across all scenarios."""
    scenarios = []

    # Timing: one row per T value at PC10
    for d in sorted(DATA_DIR.glob("msprime_timing_T*")):
        T = int(d.name.split("T")[-1])
        lams = load_lambdas(d)
        if lams:
            scenarios.append((f"Timing T={T}", "timing", lams.get(10, None)))

    for key, dirname, label in [
        ("ghost",    "msprime_ghost",               "Ghost (C absent)"),
        ("control",  "msprime_ghost_control",        "Ghost (C present)"),
        ("ss_full",  "msprime_stepping_stone_full",  "IBD full gradient"),
        ("ss_ends",  "msprime_stepping_stone_endpoints", "IBD endpoints"),
        ("pulse",    "msprime_flow_pulse",            "Pulse admixture"),
        ("cont",     "msprime_flow_continuous",       "Continuous flow"),
    ]:
        lams = load_lambdas(DATA_DIR / dirname)
        if lams:
            scenarios.append((label, key, lams.get(10, None)))

    if not scenarios:
        print("  summary: no results yet")
        return

    labels = [s[0] for s in scenarios]
    values = [s[2] for s in scenarios]
    colors_bar = ["#e41a1c" if v and v > 1.05 else "#377eb8" for v in values]

    fig = go.Figure(go.Bar(
        x=labels, y=values,
        marker_color=colors_bar,
        text=[f"{v:.3f}" if v else "N/A" for v in values],
        textposition="outside",
    ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray")
    fig.add_hline(y=1.05, line_dash="dot", line_color="orange",
                  annotation_text="λ=1.05 (notable inflation)")
    fig.update_layout(
        title="Summary: λ at PC10 Across All msprime Experiments",
        yaxis_title="Genomic inflation factor λ (PC10)",
        xaxis_tickangle=-35,
        template="plotly_white",
        width=1000, height=550,
        yaxis=dict(range=[0.8, max(v for v in values if v) * 1.15]),
    )
    out = OUT_DIR / "summary_pc10.html"
    fig.write_html(out)
    print(f"  [summary] saved {out}")


def _print_table(title: str, data: dict, row_label: str = "Scenario"):
    if not data:
        return
    all_pcs = sorted({k for lams in data.values() for k in lams if k > 0})
    header = f"\n{'─'*60}\n{title}\n{'─'*60}"
    print(header)
    col_w = 8
    print(f"{'':25}" + "".join(f"{'PC'+str(n):>{col_w}}" for n in all_pcs))
    for name, lams in data.items():
        row = f"{str(name)[:25]:<25}"
        for n in all_pcs:
            v = lams.get(n)
            row += f"{v:>{col_w}.4f}" if v is not None else f"{'—':>{col_w}}"
        print(row)


def plot_ghost_c_sweep():
    summary = DATA_DIR.parent / "outputs" / "msprime_ghost_c_sweep_summary.json"
    if not summary.exists():
        print("  ghost_c_sweep: no results yet")
        return
    with open(summary) as f:
        raw = json.load(f)
    # Keys: int PC counts + "bolt" string
    data = {}
    for k, lams in raw.items():
        parsed = {}
        for n, v in lams.items():
            try:
                parsed[int(n)] = v
            except ValueError:
                parsed[n] = v   # "bolt"
        data[int(k)] = parsed

    c_vals = sorted(data)
    fig = go.Figure()
    for i, n in enumerate(PC_COUNTS):
        y = [data[c].get(n) for c in c_vals]
        if any(v is not None for v in y):
            fig.add_trace(go.Scatter(
                x=c_vals, y=y, mode="lines+markers",
                name=f"PC{n}", line=dict(color=COLORS[i]),
                marker=dict(size=8),
            ))
    # BOLT-LMM trace
    bolt_y = [data[c].get("bolt") for c in c_vals]
    if any(v is not None for v in bolt_y):
        fig.add_trace(go.Scatter(
            x=c_vals, y=bolt_y, mode="lines+markers",
            name="BOLT-LMM", line=dict(color="black", dash="dot"),
            marker=dict(size=8, symbol="diamond"),
        ))
    fig.add_hline(y=1.0, line_dash="dash", line_color="gray", annotation_text="λ=1 (null)")
    fig.update_layout(
        title="Ghost C-Sample Sweep: λ vs Number of Pure C Samples Added<br>"
              "<sup>Base cohort: A(400) + ADMIX_AC(300) fixed. Phenotype: A=0, C=1, ADMIX_AC=1.</sup>",
        xaxis_title="Pure C samples added",
        yaxis_title="Genomic inflation factor λ",
        legend_title="Method",
        template="plotly_white",
        width=800, height=500,
    )
    out = OUT_DIR / "ghost_c_sweep.html"
    fig.write_html(out)
    print(f"  [ghost_c_sweep] saved {out}")

    # Print table
    pc_counts = sorted({k for v in data.values() for k in v if isinstance(k, int) and k > 0})
    print(f"\n{'─'*60}\nGhost C-Sample Sweep\n{'─'*60}")
    print(f"{'C samples':>10}" + "".join(f"  PC{n:>2}" for n in pc_counts) + "    BOLT")
    for c in c_vals:
        lams = data[c]
        row = f"{c:>10}" + "".join(
            f"  {lams.get(n, float('nan')):>5.3f}" for n in pc_counts
        )
        bolt_v = lams.get("bolt")
        row += f"  {bolt_v:.3f}" if bolt_v else "      —"
        print(row)


if __name__ == "__main__":
    print("=== msprime PCA Failure Results ===\n")
    plot_timing()
    plot_ghost()
    plot_ghost_c_sweep()
    plot_stepping()
    plot_flow()
    plot_summary()
    print(f"\nAll plots → {OUT_DIR}")
