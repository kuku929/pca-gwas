from pathlib import Path
from os import PathLike

from pydantic_yaml import to_yaml_file
from pydantic import BaseModel
import yaml
from typing import List, Dict, Optional, Union
from enum import StrEnum

CODE_DIR = Path(__file__).parent
DATA_DIR = CODE_DIR.parent / "data"

class Population(StrEnum):
	African = "AFR"
	American = "AMR"
	EastAsian = "EAS"
	European = "EUR"
	CentralAsian = "CSA"
	MidEastern = "MID"
	none 		= "none"

####################################
# GLOBAL PARAMETERS
####################################

class GlobalParameters(BaseModel):
	random_seed: int
	chromosome: Union[int, str]
	superpopulation: Population 
	memory: int
	batchsize: int


####################################
# FILEPATHS
####################################

class GeneralPaths(BaseModel):
	output_dir: str
	output_prefix: str


class GenotypePaths(BaseModel):
	vcf_input_raw: str
	vcf_input_processed: str
	vcf_metadata: str
	popfile_raw: str
	popfile_processed: str
	variant_list: str
	remove_list: str
	rsid_list: str
	genetic_mapfile: str
	genetic_distfile: str
	mutation_mapfile: str
	mutation_agefile: str
	hap1_matrix: str
	hap2_matrix: str


class PhenotypePaths(BaseModel):
	causal_list: str
	reference: str
	plink_override: Optional[str]


class SoftwarePaths(BaseModel):
	plink: str
	plink2: str
	king: str
	vcftools: str
	mapthin: str
	phenoalg: str


class Filepaths(BaseModel):
	general: GeneralPaths
	genotype: GenotypePaths
	phenotype: PhenotypePaths
	software: SoftwarePaths


####################################
# GENOTYPE DATA
####################################

class PopulationGroup(BaseModel):
	id: str
	nsamples: int
	populations: List[Dict[str, int]]


class Samples(BaseModel):
	use_default: bool
	custom: List[PopulationGroup]
	default: Dict[str, int]


class GenotypeData(BaseModel):
	samples: Samples
	rho: Dict[str, float]
	Ne: Dict[str, int]


####################################
# PHENOTYPE DATA
####################################

class Causality(BaseModel):
	UseCausalList: bool
	Polygenicity: Optional[float]
	Pleiotropy: Optional[float]


class PhenotypeData(BaseModel):
	nPopulation: int
	nTrait: int
	a: float
	b: float
	c: float
	nComponent: int
	PropotionGeno: str
	PropotionCovar: str
	Prevalence: str
	TraitCorr: float
	PopulationCorr: str
	CompWeight: str
	Causality: Causality


####################################
# EVALUATION
####################################

class Metrics(BaseModel):
	aats: bool
	kinship: bool
	ld_corr: bool
	ld_decay: bool
	maf: bool
	pca: bool
	gwas: bool


class Evaluation(BaseModel):
	metrics: Metrics


####################################
# OPTIMISATION
####################################

class UniformPrior(BaseModel):
	uniform_lower: float
	uniform_upper: float


class Priors(BaseModel):
	rho: UniformPrior
	Ne: UniformPrior


class SimulationABC(BaseModel):
	run: bool
	n_particles: int
	threshold: float
	max_iter: int
	write_progress: bool


class EmulationABC(BaseModel):
	run: bool
	n_particles: int
	threshold: float
	n_design_points: int
	max_iter: int
	write_progress: bool


class SummaryStatistics(BaseModel):
	ld_decay: bool
	kinship: bool


class Optimisation(BaseModel):
	priors: Priors
	simulation_rejection_ABC: SimulationABC
	emulation_rejection_ABC: EmulationABC
	summary_statistics: SummaryStatistics


####################################
# ROOT CONFIG
####################################

class Config(BaseModel):
	global_parameters: GlobalParameters
	filepaths: Filepaths
	genotype_data: GenotypeData
	phenotype_data: PhenotypeData
	evaluation: Evaluation
	optimisation: Optimisation

	def default(self):
		self.write()

	def homogenous_population(self, pop : Population, nsamples : int):
		self.global_parameters.superpopulation = pop
		self.genotype_data.samples.use_default = True
		self.genotype_data.samples.default = {"nsamples" : nsamples}

	def mixed_population(self, groups : List[PopulationGroup]):
		self.global_parameters.superpopulation = Population.none
		self.genotype_data.samples.use_default = False
		self.genotype_data.samples.custom = groups

	def set_phenotype_params(self, n_populations: int, heritability: Union[float, List[float]] = 0.1):
		identity_row = lambda i, n: ",".join("1" if i == j else "0" for j in range(n))
		pop_corr = ",".join(identity_row(i, n_populations) for i in range(n_populations))
		if isinstance(heritability, (int, float)):
			h2_list = [heritability] * n_populations
		else:
			h2_list = list(heritability)
		self.phenotype_data.nPopulation    = n_populations
		self.phenotype_data.PropotionGeno  = ",".join(str(h) for h in h2_list)
		self.phenotype_data.PropotionCovar = ",".join(["0"]   * n_populations)
		self.phenotype_data.Prevalence     = ",".join(["0.5"] * n_populations)
		self.phenotype_data.PopulationCorr = pop_corr
	
	def set_output_dir(self, suffix : str):
		self.filepaths.general.output_dir = "data/outputs/test_" + suffix
		print("=> Output dir: ", self.filepaths.general.output_dir)
	
	def _write_to(self, path : PathLike):
		print("=> Writing to: ", path)
		try:
			to_yaml_file(path / "config.yaml", self)
		except FileNotFoundError as e:
			print(e)
			print("Did you generate the genotypes?")
		except PermissionError as e:
			print(f"Warning: permission denied writing config to {path}: {e}")

	def write(self):
		self._write_to(DATA_DIR)
	
	def save_config(self):
		self._write_to(DATA_DIR.parent / self.filepaths.general.output_dir)

with open(CODE_DIR / "template.yaml", "r") as fp:
	data = yaml.safe_load(fp)

config = Config(**data)
