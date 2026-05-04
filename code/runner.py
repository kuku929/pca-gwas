"""
Runs HAPNEST simulation with the generated config.
"""
from pathlib import Path
from os import PathLike
import io
import re
import shutil
import sys
import subprocess
import time
import json
import numpy as np
import pandas as pd
from scipy import stats as scipy_stats
from configure import config
CODE_DIR = Path(__file__).parent
HOST_DATA_DIR = CODE_DIR.parent / "data"   # resolves through the symlink on the host
DOCKER_DATA_DIR = Path("data") # data dir location inside the docker container
DOCKER = ["docker", "exec", "--user", "1000:1000", "-e", "JULIA_DEPOT_PATH=/root/.julia", "hapnest"]

def run(command : str, file : PathLike):
	try:
		with io.open(file, "wb", buffering=0) as writer, io.open(file, "rb", buffering=0) as reader:
			process = subprocess.Popen(command, stdout=writer, stderr=subprocess.STDOUT, preexec_fn=lambda : print(command, flush=True))
			while process.poll() is None:
				sys.stdout.buffer.write(reader.read())
				time.sleep(0.5)
				# Read the remaining
				sys.stdout.buffer.write(reader.read())
	except Exception as e:
		print(e)
		print("Did you start the docker container? Use this command(make sure your path is correct while executing this):" \
		"NOTE: the docker container must be named 'hapnest'.")

def generate_geno(nthreads=1):
	command = DOCKER + ["generate_geno", str(nthreads), (DOCKER_DATA_DIR / "config.yaml").as_posix()]
	file = CODE_DIR / "log/geno.out"
	run(command, file)
	print("=> Done generating genotype")

def generate_pheno():
	command = DOCKER + ["generate_pheno", (DOCKER_DATA_DIR / "config.yaml").as_posix()]
	file = CODE_DIR / "log/pheno.out"
	run(command, file)
	print("=> Done generating phenotype")

def run_gwas_pc_sweep(output_dir: str, prefix_template: str, chromosome, pc_counts=[0, 2, 5, 10]):
	"""
	Run plink2 GWAS inside Docker with 0, 2, 5, 10 PCs as covariates.
	Reuses existing experiment data — no generate_geno/pheno needed.
	Saves lambda values to evaluation/pc_sweep_lambdas.json and returns them.
	"""
	base       = re.sub(r'-?\{chromosome\}', '', prefix_template)   # "test_chr"
	chr_pfx    = f"{base}-{chromosome}"                             # "test_chr-2"
	host_out   = HOST_DATA_DIR.parent / output_dir
	eval_dir   = host_out / "evaluation"
	docker_out = Path(output_dir)

	# Build temp phenotype file (FID IID PHENO) from .pheno1
	# generate_pheno writes {base}.pheno1 (no chromosome suffix)
	pheno_src = host_out / f"{base}.pheno1"
	pheno_df  = pd.read_csv(pheno_src, sep='\t')
	sample0   = str(pheno_df['Sample'].iloc[0])
	if '_' in sample0:
		fid = pheno_df['Sample'].str.rsplit('_', n=1).str[0]
		iid = pheno_df['Sample'].str.rsplit('_', n=1).str[1]
	else:
		fid = iid = pheno_df['Sample']
	liability_col = next(c for c in pheno_df.columns if 'liability' in c.lower())
	pheno_tmp = eval_dir / "pheno_sweep.tmp"
	pd.DataFrame({'#FID': fid, 'IID': iid, 'PHENO': pheno_df[liability_col]}).to_csv(
		pheno_tmp, sep='\t', index=False
	)

	# Load PCA eigenvec for covariate files
	eigenvec = eval_dir / f"{chr_pfx}.syn.pca.eigenvec"
	pca_df   = pd.read_csv(eigenvec, sep=r'\s+')
	pca_df.columns = ['FID', 'IID'] + [f'PC{i}' for i in range(1, len(pca_df.columns) - 1)]

	lambdas = {}
	for n in pc_counts:
		tag        = f"gwas_pc{n}"
		out_docker = docker_out / "evaluation" / f"{chr_pfx}.{tag}"

		cmd = DOCKER + [
			"plink2",
			"--bed",        (docker_out / f"{chr_pfx}.bed").as_posix(),
			"--bim",        (docker_out / f"{chr_pfx}.bim").as_posix(),
			"--fam",        (docker_out / f"{chr_pfx}.fam").as_posix(),
			"--pheno",      (docker_out / "evaluation" / "pheno_sweep.tmp").as_posix(),
			"--pheno-name", "PHENO",
			"--glm",        "hide-covar",
			"--maf",        "0.01",
			"--allow-no-sex",
			"--out",        out_docker.as_posix(),
		]
		if n > 0:
			covar_host = eval_dir / f"covar_pc{n}.tmp"
			covar_df   = pca_df[['FID', 'IID'] + [f'PC{i}' for i in range(1, n + 1)]].copy()
			covar_df   = covar_df.rename(columns={'FID': '#FID'})
			covar_df.to_csv(covar_host, sep='\t', index=False)
			cmd += ["--covar", (docker_out / "evaluation" / f"covar_pc{n}.tmp").as_posix()]

		run(cmd, CODE_DIR / "log" / f"{tag}.out")

		glm_files = sorted((eval_dir).glob(f"{chr_pfx}.{tag}*.glm.linear"))
		if not glm_files:
			print(f"  PC{n}: no output found — check log/{tag}.out")
			continue
		glm = pd.read_csv(glm_files[0], sep='\t')
		glm.columns = [c.lstrip('#') for c in glm.columns]
		if 'TEST' in glm.columns:
			glm = glm[glm['TEST'] == 'ADD']
		pvals = pd.to_numeric(glm['P'], errors='coerce').dropna().values
		pvals = pvals[pvals > 0]
		lam   = float(np.median(scipy_stats.chi2.ppf(1 - pvals, df=1)) / 0.456)
		lambdas[n] = round(lam, 4)
		print(f"  PC{n:2d}: λ = {lam:.4f}")

	with open(eval_dir / "pc_sweep_lambdas.json", 'w') as f:
		json.dump(lambdas, f, indent=2)
	print(f"=> saved {eval_dir / 'pc_sweep_lambdas.json'}")
	return lambdas

def validate():
	command = DOCKER + ["validate", (DOCKER_DATA_DIR / "config.yaml").as_posix()]
	file = CODE_DIR / "log/validate.out"
	run(command, file)
	try:
		run_gwas_pc_sweep(
				output_dir=config.filepaths.general.output_dir,
				prefix_template="test_chr-{chromosome}",
				chromosome=2,
				pc_counts=[0, 2, 5, 10],
		)
	except Exception as e:
		print(e)
		print("Run: docker exec hapnest chmod 777 -R /data")
	print("=> Validation complete")
