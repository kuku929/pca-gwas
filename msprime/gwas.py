"""plink2 GWAS pipeline (via Docker) for msprime-generated data."""
import json
import subprocess
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from pathlib import Path

PROJECT_ROOT = Path(__file__).parent.parent
DATA_DIR = PROJECT_ROOT / "data"
DOCKER_DATA = Path("data")
DOCKER = ["docker", "exec", "--user", "1000:1000", "hapnest"]
BOLT_BIN = PROJECT_ROOT / "BOLT-LMM_v2.5" / "bolt"


def _to_docker(host_path: Path) -> Path:
    return DOCKER_DATA / host_path.relative_to(DATA_DIR)


def _run(cmd: list, label: str = "plink2"):
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0:
        print(f"[{label}] failed (rc={result.returncode})")
        if result.stderr:
            print(result.stderr[-800:])
    return result.returncode == 0


def vcf_to_bed(vcf_host: Path, bed_prefix_host: Path, delete_vcf: bool = True) -> bool:
    bed_prefix_host.parent.mkdir(parents=True, exist_ok=True)
    ok = _run(DOCKER + [
        "plink2",
        "--vcf", str(_to_docker(vcf_host)),
        "--double-id",
        "--max-alleles", "2",
        "--maf", "0.01",
        "--make-bed",
        "--allow-no-sex",
        "--out", str(_to_docker(bed_prefix_host)),
    ], "vcf_to_bed")
    if ok:
        print(f"BED: {bed_prefix_host}.bed")
        if delete_vcf and vcf_host.exists():
            vcf_host.unlink()
            print(f"VCF deleted (space saved)")
    return ok


def run_pca(bed_prefix_host: Path, n_pcs: int = 20) -> bool:
    ok = _run(DOCKER + [
        "plink2",
        "--bfile", str(_to_docker(bed_prefix_host)),
        "--pca", str(n_pcs),
        "--allow-no-sex",
        "--out", str(_to_docker(bed_prefix_host)),
    ], "pca")
    if ok:
        print(f"PCA: {bed_prefix_host}.eigenvec")
    return ok


def gwas_pc_sweep(
    bed_prefix_host: Path,
    pheno_host: Path,
    eigenvec_host: Path,
    pc_counts: list = None,
) -> dict:
    """Run GWAS with each PC count. Returns {n_pc: lambda}."""
    if pc_counts is None:
        pc_counts = [0, 2, 5, 10, 20]

    pca_df = pd.read_csv(eigenvec_host, sep=r"\s+")
    pca_df.columns = ["FID", "IID"] + [f"PC{i}" for i in range(1, len(pca_df.columns) - 1)]

    lambdas = {}
    out_dir = bed_prefix_host.parent
    for n in pc_counts:
        out_prefix = out_dir / f"gwas_pc{n}"
        cmd = DOCKER + [
            "plink2",
            "--bfile",      str(_to_docker(bed_prefix_host)),
            "--pheno",      str(_to_docker(pheno_host)),
            "--pheno-name", "PHENO",
            "--glm",        "hide-covar", "allow-no-covars",
            "--maf",        "0.01",
            "--allow-no-sex",
            "--out",        str(_to_docker(out_prefix)),
        ]
        if n > 0:
            covar_path = out_dir / f"covar_pc{n}.txt"
            covar_df = pca_df[["FID", "IID"] + [f"PC{i}" for i in range(1, n + 1)]].copy()
            covar_df = covar_df.rename(columns={"FID": "#FID"})
            covar_df.to_csv(covar_path, sep="\t", index=False)
            cmd += ["--covar", str(_to_docker(covar_path))]

        ok = _run(cmd, f"gwas_pc{n}")
        if not ok:
            continue

        glm_files = sorted(out_dir.glob(f"gwas_pc{n}*.glm.linear"))
        if not glm_files:
            print(f"  PC{n}: no .glm.linear output")
            continue
        glm = pd.read_csv(glm_files[0], sep="\t")
        glm.columns = [c.lstrip("#") for c in glm.columns]
        if "TEST" in glm.columns:
            glm = glm[glm["TEST"] == "ADD"]
        pvals = pd.to_numeric(glm["P"], errors="coerce").dropna().values
        pvals = pvals[pvals > 0]
        if len(pvals) == 0:
            print(f"  PC{n}: no valid p-values")
            continue
        lam = float(np.median(scipy_stats.chi2.ppf(1 - pvals, df=1)) / 0.456)
        lambdas[n] = round(lam, 4)
        print(f"  PC{n:2d}: λ = {lam:.4f}")

    results_path = out_dir / "lambdas.json"
    with open(results_path, "w") as f:
        json.dump(lambdas, f, indent=2)
    print(f"=> saved {results_path}")
    return lambdas


def run_bolt(bed_prefix_host: Path, pheno_host: Path) -> float | None:
    """Run BOLT-LMM. Returns λ or None on failure."""
    if not BOLT_BIN.exists():
        print(f"  BOLT-LMM binary not found: {BOLT_BIN}")
        return None
    stats_file = bed_prefix_host.parent / "bolt.stats"
    cmd = [
        str(BOLT_BIN),
        "--bed",               str(bed_prefix_host) + ".bed",
        "--bim",               str(bed_prefix_host) + ".bim",
        "--fam",               str(bed_prefix_host) + ".fam",
        "--phenoFile",         str(pheno_host),
        "--phenoCol",          "PHENO",
        "--lmmForceNonInf",
        "--LDscoresUseChip",
        "--numLeaveOutChunks", "2",
        "--statsFile",         str(stats_file),
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    if result.returncode != 0 or not stats_file.exists():
        print(f"  BOLT-LMM failed (rc={result.returncode})")
        if result.stderr:
            print(result.stderr[-600:])
        return None
    df_stats = pd.read_csv(stats_file, sep="\t")
    p_col = "P_BOLT_LMM" if "P_BOLT_LMM" in df_stats.columns else "P_BOLT_LMM_INF"
    pvals = pd.to_numeric(df_stats[p_col], errors="coerce").dropna().values
    pvals = pvals[pvals > 0]
    lam = float(np.median(scipy_stats.chi2.ppf(1 - pvals, df=1)) / 0.456)
    print(f"  BOLT-LMM: λ = {lam:.4f}")
    return round(lam, 4)


def run_experiment(
    output_dir: Path,
    individual_names: list,
    pop_labels: list,
    phenotype: np.ndarray,
    pc_counts: list = None,
    run_bolt_lmm: bool = False,
) -> dict:
    """Full pipeline: VCF → BED → PCA → GWAS sweep. Returns lambda dict."""
    vcf_host = output_dir / "sim.vcf"
    bed_prefix = output_dir / "sim"
    eigenvec_host = output_dir / "sim.eigenvec"

    from phenotype import write_pheno
    pheno_host = write_pheno(output_dir, individual_names, phenotype)

    if not vcf_to_bed(vcf_host, bed_prefix):
        raise RuntimeError("vcf_to_bed failed — is hapnest Docker running?")
    if not run_pca(bed_prefix, n_pcs=max(pc_counts or [20])):
        raise RuntimeError("PCA failed")

    lambdas = gwas_pc_sweep(bed_prefix, pheno_host, eigenvec_host, pc_counts)

    if run_bolt_lmm:
        # BOLT-LMM needs FID col without # — write a separate pheno file
        bolt_pheno = output_dir / "bolt_pheno.txt"
        df_p = pd.read_csv(pheno_host, sep="\t")
        df_p.columns = [c.lstrip("#") for c in df_p.columns]
        df_p.to_csv(bolt_pheno, sep="\t", index=False)
        lam_bolt = run_bolt(bed_prefix, bolt_pheno)
        if lam_bolt is not None:
            lambdas["bolt"] = lam_bolt
            with open(output_dir / "lambdas.json", "w") as f:
                json.dump(lambdas, f, indent=2)

    return lambdas
