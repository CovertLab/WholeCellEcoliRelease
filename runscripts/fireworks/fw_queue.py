#!/usr/bin/env python

'''
Creates an array of firetasks (wf_fws) and specifies their links (wf_links) in
a workflow (wf) for Fireworks, and submits them to the queue.

Several environmental variables can be specified, shown below with their (type, default value).
These are set as follows, and otherwise revert to their default value:

	VARIABLE=VALUE python runscripts/fireworks/fw_queue.py

Set description:
	DESC (str, ""): a description of the simulation, used to name output folder

Variant variables:
	VARIANT (str, "wildtype"): specifies the environmental condition, see
		models/ecoli/sim/variants/*.py for the variant choices
	FIRST_VARIANT_INDEX (int, "0"): index of the first variant condition to run;
		the control index depends on the particular variant condition
	LAST_VARIANT_INDEX (int, "0"): index of the last variant condition to run; Fireworks
		will run all conditions between FIRST_VARIANT_INDEX and LAST_VARIANT_INDEX

Workflow options:
	CACHED_SIM_DATA (int, "0"): if nonzero, previously cached data will be used
		to run the simulation instead of running the parca; useful for repeated
		simulations where raw_data/parca files are not changed
	PARALLEL_PARCA (int, "0"): if nonzero, some parca operations will run in
		parallel instead of serially
	DEBUG_PARCA (int, "0"): if nonzero, this makes Parca calculate only one
		arbitrarily-chosen transcription factor condition when adjusting gene
		expression levels, leaving the others at their input levels, for faster
		Parca debugging; do not use this for an actual simulation
	COMPRESS_OUTPUT (int, "0"): if nonzero, outputs will be compressed (.bz2)
	RUN_AGGREGATE_ANALYSIS (int, "1"): if nonzero, all analyses (aggregate and
		single gen) are run on simulation output
	PLOTS (str, ""): Which analyses to run (if RUN_AGGREGATE_ANALYSIS is
		true). This should name one or more tags. For more than one tag,
		separate them with whitespace and remember shell quoting. ACTIVE
		includes all active plots. CORE includes just the plots recommended for
		everyday development. VARIANT runs analysis specific to the specified
		variant. "DEFAULT" runs both "CORE" and "VARIANT" and is selected by default.
		You can also name specific analysis files but any analysis
		categories that don't have such a filename will print error messages.
	DISABLE_RIBOSOME_CAPACITY_FITTING (int, "0"): if nonzero, ribosome
		expression is not fit to protein synthesis demands
	DISABLE_RNAPOLY_CAPACITY_FITTING (int, "0"): if nonzero, RNA polymerase
		expression is not fit to RNA synthesis demands
	WC_ANALYZE_FAST (anything, --): if set, run each analysis plot in a separate
		process
	BUILD_CAUSALITY_NETWORK (int, "0"): if nonzero, causality network files are
		generated from the first sim's output
	RAISE_ON_TIME_LIMIT (int, "0"), if nonzero, the simulation raises an error
		if the time limit (WC_LENGTHSEC) is reached before division

Simulation parameters:
	N_GENS (int, "1"): the number of generations to be simulated
	N_INIT_SIMS (int, "1"): the number of initial simulations
	SEED (int, "0"): starting seed to run
	SINGLE_DAUGHTERS (int, "1"): if nonzero, the simulation will generate only
		one daughter cell for each new generation rather than two, thus avoiding
		an exponential increase in the number of simulations
	TIMELINE (str, "0 minimal"): sets the timeline of events for the simulation.
		See wholecell/utils/make_media.py, make_timeline() for timeline
		formatting details.
	WC_LENGTHSEC (int, "10800"): sets the maximum simulation time in seconds, useful
		for short simulations (default is 3 hr)
	TIMESTEP_MAX (float, "2"): sets the maximum time step
	TIMESTEP_SAFETY_FRAC (float, "1.3"): increases the time step by this factor
		if conditions are favorable; up the the limit of the max time step
	TIMESTEP_UPDATE_FREQ (int, "5"): frequency at which the time step is updated
	ADJUST_TIMESTEP_FOR_CHARGING (int, "0"): if nonzero, adjusts the timestep
		if charging creates a large update to improve stability of sims
	LOG_TO_DISK_EVERY (int, "1"): frequency at which simulation outputs are
		logged to disk
	JIT (int, "1"): if nonzero, jit compiled functions are used for certain
		processes, otherwise only uses lambda functions

Modeling options:
	MASS_DISTRIBUTION (int, "1"): if nonzero, a mass coefficient is drawn from
		a normal distribution centered on 1; otherwise it is set equal to 1
	GROWTH_RATE_NOISE (int, "0"): if nonzero, a growth rate coefficient is drawn
		from a normal distribution centered on 1; otherwise it is set equal to 1
	D_PERIOD_DIVISION (int, "0"): if nonzero, ends simulation once D period has
		occurred after chromosome termination; otherwise simulation terminates
		once a given mass has been added to the cell
	OPERONS (str, "off"): run with operons "off", "on" (actually monocistronic
		or polycistronic), or "both" (into adjacent output directories)
	VARIABLE_ELONGATION_TRANSCRIPTION (int, "1"): if nonzero, use variable
		transcription elongation rates for each gene
	VARIABLE_ELONGATION_TRANSLATION (int, "0"): if nonzero, use variable
		translation elongation rates for each gene
	TRANSLATION_SUPPLY (int, "1"): if nonzero, the ribosome elongation rate is
		limited by the condition specific rate of amino acid supply; otherwise
		the elongation rate is set by condition
	TRNA_CHARGING (int, "1"): if nonzero, tRNA charging reactions are modeled
		and the ribosome elongation rate is set by the amount of charged tRNA
		present.  This option will override TRANSLATION_SUPPLY in the simulation.
	AA_SUPPLY_IN_CHARGING (int, "0"): if nonzero, amino acid supply function is
		used during charging for more stable charging calculations.  Only has an
		effect if TRNA_CHARGING option is used.
	PPGPP_REGULATION (int, "0"): if nonzero, ppGpp concentration is determined
		with kinetic equations
	SUPERHELICAL_DENSITY (int, "0"): if nonzero, dynamically compute
		superhelical densities of each DNA fragment
	RECYCLE_STALLED_ELONGATION (int "0"): if nonzero, recycle RNAP and fragment
		bases when transcription elongation is stalled in ntp-limiting conditions
	MECHANISTIC_REPLISOME (int, "1"): if nonzero, replisome initiation is
		mechanistic (requires appropriate number of subunits to initiate)
	MECHANISTIC_TRANSLATION_SUPPLY (int, "0"): if nonzero, amino acid translation supply is
		mechanistic (depends on concentrations of enzymes and amino acids)
	MECHANISTIC_AA_TRANSPORT (int, "0"): if nonzero, amino acid uptake and export are
		mechanistic (depends on concentrations of transporter enzymes and amino acids)
	TRNA_ATTENUATION (int, "0"): if nonzero, transcriptional attenuation by
		charged tRNA is enabled

Additional variables:
	LAUNCHPAD_FILE (str, "my_launchpad.yaml"): set launchpad config file location
	VERBOSE_QUEUE (int, "1"): if nonzero, gives more detailed messages during
		fireworks set up

Environment variables that matter when running the workflow:
	DEBUG_GC (int, "0"): if nonzero, enable leak detection in the analysis plots


-------------------------------------------------
Current dependency network of existing Firetasks (each task has a list of its
immediate downstream dependencies):

InitRawData (fw_init_raw_data)
* FitSimData
* InitValidationData

InitRawValidationData (fw_raw_validation_data)
* InitValidationData

InitValidationData (fw_validation_data)
* CompressRawValidationData if COMPRESS_OUTPUT
* CompressRawData if COMPRESS_OUTPUT
* AnalysisParca if RUN_AGGREGATE_ANALYSIS
* VariantSimData [see the comment in the code]

FitSimData (fw_calculate_sim_data)
* CompressRawData if COMPRESS_OUTPUT
* VariantSimData * VARIANTS_TO_RUN
* AnalysisParca if RUN_AGGREGATE_ANALYSIS
* CompressFitSimData if COMPRESS_OUTPUT

VariantSimData (fw_this_variant_sim_data)
* CompressFitSimData if COMPRESS_OUTPUT
* Simulation * N_INIT_SIMS

Simulation/SimulationDaughter (fw_this_variant_this_gen_this_sim)
* SimulationDaughter (* 2 if not SINGLE_DAUGHTERS) if GEN < N_GENS - 1
* AnalysisSingle if RUN_AGGREGATE_ANALYSIS
* AnalysisCohort if RUN_AGGREGATE_ANALYSIS and GEN == N_GENS - 1
* AnalysisVariant if RUN_AGGREGATE_ANALYSIS and GEN == N_GENS - 1
* AnalysisMultiGen if RUN_AGGREGATE_ANALYSIS and GEN == N_GENS - 1
* BuildCausalityNetwork if BUILD_CAUSALITY_NETWORK

AnalysisParca (fw_parca_analysis)
* CompressFitSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT

AnalysisSingle (fw_this_variant_this_gen_this_sim_analysis)
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

AnalysisVariant (fw_variant_analysis)
* AnalysisComparison if OPERONS == 'both'
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

AnalysisCohort (fw_this_variant_cohort_analysis)
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

AnalysisMultiGen (fw_this_variant_this_seed_multigen_analysis)
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

BuildCausalityNetwork (fw_this_variant_this_gen_this_sim_causality_network)
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

CompressValidationData (fw_validation_data_compression)

CompressRawData (fw_raw_data_compression)

CompressFitSimData (fw_sim_data_1_compression)

CompressRawValidationData (fw_raw_validation_data_compression)

CompressVariantSimData (fw_this_variant_sim_data_compression)

CompressSimulationOutput (fw_this_variant_this_gen_this_sim_compression)
-------------------------------------------------

'''

from __future__ import annotations

import collections
import os
import sys
from typing import Any, Dict, List, Optional, Union

import yaml
from fireworks import FiretaskBase, Firework, LaunchPad, Workflow, ScriptTask

from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.fireworks.firetasks import SimulationTask
from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.fireworks.firetasks import AnalysisParcaTask
from wholecell.fireworks.firetasks import AnalysisVariantTask
from wholecell.fireworks.firetasks import AnalysisCohortTask
from wholecell.fireworks.firetasks import AnalysisSingleTask
from wholecell.fireworks.firetasks import AnalysisMultiGenTask
from wholecell.fireworks.firetasks import AnalysisComparisonTask
from wholecell.fireworks.firetasks import BuildCausalityNetworkTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils import constants
from wholecell.utils import filepath


def get_environment(variable, default):
	'''
	Gets an environmental variable value and prints a confirmation of the
	variable and value if it is in the environment.

	Args:
		variable (str): name of environmental variable to get
		default (str): default value if variable is not in the environment

	Returns:
		str: value of the environmental variable or default if it doesn't exist
	'''

	value = os.environ.get(variable, None)

	if value is not None:
		print('\t{}: {}'.format(variable, value))
	else:
		value = default

	return value

### Read input options

print('Parsed environmental variables:')

### Set variant variables

VARIANT = get_environment("VARIANT", "wildtype")
FIRST_VARIANT_INDEX = int(get_environment("FIRST_VARIANT_INDEX", "0"))
LAST_VARIANT_INDEX = int(get_environment("LAST_VARIANT_INDEX", "0"))

# This variable gets iterated over in multiple places
# So be careful if you change it to xrange
VARIANTS_TO_RUN = list(range(FIRST_VARIANT_INDEX, LAST_VARIANT_INDEX + 1))

### Set other simulation parameters

TIMELINE = str(get_environment("TIMELINE", DEFAULT_SIMULATION_KWARGS["timeline"]))
WC_LENGTHSEC = int(get_environment("WC_LENGTHSEC", DEFAULT_SIMULATION_KWARGS["lengthSec"]))
TIMESTEP_SAFETY_FRAC = float(get_environment("TIMESTEP_SAFETY_FRAC", DEFAULT_SIMULATION_KWARGS["timeStepSafetyFraction"]))
TIMESTEP_MAX = float(get_environment("TIMESTEP_MAX", DEFAULT_SIMULATION_KWARGS["maxTimeStep"]))
TIMESTEP_UPDATE_FREQ = int(get_environment("TIMESTEP_UPDATE_FREQ", DEFAULT_SIMULATION_KWARGS["updateTimeStepFreq"]))
ADJUST_TIMESTEP_FOR_CHARGING = int(get_environment("ADJUST_TIMESTEP_FOR_CHARGING", DEFAULT_SIMULATION_KWARGS["adjust_timestep_for_charging"]))
LOG_TO_DISK_EVERY = int(get_environment("LOG_TO_DISK_EVERY", DEFAULT_SIMULATION_KWARGS["logToDiskEvery"]))
JIT = bool(int(get_environment("JIT", DEFAULT_SIMULATION_KWARGS["jit"])))
MASS_DISTRIBUTION = bool(int(get_environment("MASS_DISTRIBUTION", DEFAULT_SIMULATION_KWARGS["massDistribution"])))
GROWTH_RATE_NOISE = bool(int(get_environment("GROWTH_RATE_NOISE", DEFAULT_SIMULATION_KWARGS["growthRateNoise"])))
D_PERIOD_DIVISION = bool(int(get_environment("D_PERIOD_DIVISION", DEFAULT_SIMULATION_KWARGS["dPeriodDivision"])))
OPERONS = get_environment("OPERONS", constants.DEFAULT_OPERON_OPTION)
assert OPERONS in constants.EXTENDED_OPERON_OPTIONS, f'{OPERONS=} needs to be in {constants.EXTENDED_OPERON_OPTIONS}'
VARIABLE_ELONGATION_TRANSCRIPTION = bool(int(get_environment("VARIABLE_ELONGATION_TRANSCRIPTION", DEFAULT_SIMULATION_KWARGS["variable_elongation_transcription"])))
VARIABLE_ELONGATION_TRANSLATION = bool(int(get_environment("VARIABLE_ELONGATION_TRANSLATION", DEFAULT_SIMULATION_KWARGS["variable_elongation_translation"])))
TRANSLATION_SUPPLY = bool(int(get_environment("TRANSLATION_SUPPLY", DEFAULT_SIMULATION_KWARGS["translationSupply"])))
TRNA_CHARGING = bool(int(get_environment("TRNA_CHARGING", DEFAULT_SIMULATION_KWARGS["trna_charging"])))
AA_SUPPLY_IN_CHARGING = bool(int(get_environment("AA_SUPPLY_IN_CHARGING", DEFAULT_SIMULATION_KWARGS["aa_supply_in_charging"])))
PPGPP_REGULATION = bool(int(get_environment("PPGPP_REGULATION", DEFAULT_SIMULATION_KWARGS["ppgpp_regulation"])))
SUPERHELICAL_DENSITY = bool(int(get_environment("SUPERHELICAL_DENSITY", DEFAULT_SIMULATION_KWARGS["superhelical_density"])))
RECYCLE_STALLED_ELONGATION = bool(int(get_environment("RECYCLE_STALLED_ELONGATION", DEFAULT_SIMULATION_KWARGS["recycle_stalled_elongation"])))
MECHANISTIC_REPLISOME = bool(int(get_environment("MECHANISTIC_REPLISOME", DEFAULT_SIMULATION_KWARGS["mechanistic_replisome"])))
MECHANISTIC_TRANSLATION_SUPPLY = bool(int(get_environment("MECHANISTIC_TRANSLATION_SUPPLY", DEFAULT_SIMULATION_KWARGS["mechanistic_translation_supply"])))
MECHANISTIC_AA_TRANSPORT = bool(int(get_environment("MECHANISTIC_AA_TRANSPORT", DEFAULT_SIMULATION_KWARGS["mechanistic_aa_transport"])))
TRNA_ATTENUATION = bool(int(get_environment("TRNA_ATTENUATION", DEFAULT_SIMULATION_KWARGS["trna_attenuation"])))
RAISE_ON_TIME_LIMIT = bool(int(get_environment("RAISE_ON_TIME_LIMIT", DEFAULT_SIMULATION_KWARGS["raise_on_time_limit"])))
N_INIT_SIMS = int(get_environment("N_INIT_SIMS", "1"))
SEED = int(get_environment("SEED", "0"))
N_GENS = int(get_environment("N_GENS", "1"))
SINGLE_DAUGHTERS = bool(int(get_environment("SINGLE_DAUGHTERS", "1")))
LAUNCHPAD_FILE = str(get_environment("LAUNCHPAD_FILE", "my_launchpad.yaml"))
COMPRESS_OUTPUT = bool(int(get_environment("COMPRESS_OUTPUT", "0")))
SIM_DESCRIPTION = get_environment("DESC", "").replace(" ", "_")
VERBOSE_QUEUE = bool(int(get_environment("VERBOSE_QUEUE", "1")))
RUN_AGGREGATE_ANALYSIS = bool(int(get_environment("RUN_AGGREGATE_ANALYSIS", "1")))
PLOTS = get_environment("PLOTS", "").split()
CACHED_SIM_DATA = bool(int(get_environment("CACHED_SIM_DATA", "0")))
PARALLEL_PARCA = bool(int(get_environment("PARALLEL_PARCA", "0")))
DEBUG_PARCA = bool(int(get_environment("DEBUG_PARCA", "0")))
DISABLE_RIBOSOME_CAPACITY_FITTING = bool(int(get_environment("DISABLE_RIBOSOME_CAPACITY_FITTING", "0")))
DISABLE_RNAPOLY_CAPACITY_FITTING = bool(int(get_environment("DISABLE_RNAPOLY_CAPACITY_FITTING", "0")))
BUILD_CAUSALITY_NETWORK = bool(int(get_environment("BUILD_CAUSALITY_NETWORK", "0")))

if not RUN_AGGREGATE_ANALYSIS:
	COMPRESS_OUTPUT = False


def log_info(message):
	if VERBOSE_QUEUE:
		print(message)


OUT_DIRECTORY = filepath.makedirs(filepath.ROOT_PATH, "out")
CACHED_SIM_DATA_DIRECTORY = os.path.join(filepath.ROOT_PATH, "cached")
SUBMISSION_TIME = filepath.timestamp()

# To run each analysis plot in a separate process, ask the analysis Firetasks
# for several CPUs (it will clip to the number available) and allocate multiple
# CPUs from SLURM via the Fireworks Queue Adaptor. Otherwise, let the qadapter
# YAML file request the number of CPUs so we can tune it to ask for extra CPUs
# in order to get proportionally more RAM, e.g. after running many generations.
if "WC_ANALYZE_FAST" in os.environ:
	analysis_cpus = 8
else:
	analysis_cpus = 1


class WorkflowBuilder:
	def __init__(self) -> None:
		self.wf_fws: List[Firework] = []  # Fireworks in this workflow
		self.wf_links: Dict[Firework, List[Firework]] = collections.defaultdict(list)  # dependencies
		self.operons = ''
		self.name_suffix = ''

		# AnalysisComparisonTask depends on AnalysisVariantTask in order to
		# depend (indirectly but light weight) on all the sim tasks.
		self.fw_variant_analysis = None

	def add_firework(self, firetask: FiretaskBase,
					 name: str, *,
					 parents: Union[None, Firework, List[Firework]] = None,
					 cpus: int = 1,
					 priority: Optional[int] = None,
					 indent: int = 0) -> Firework:
		"""Construct a Firework and add it to the accumulating list.

		Args:
			firetask: add a Firework containing this Firetask.
			name: the unique task name. (In the operons=both case, this gets a
				suffix to distinguish the 'on' case, although it might not need
				to be unique.)
			parents: parent Firework(s), that is, dependencies.
			cpus: the number of CPUs to allocate to the job.
			priority: the job priority; a larger number is higher priority; None
				is lowest priority (-∞).
			indent: indentation level for the "Queueing" log message.
		"""
		TAB = '\t'
		name += self.name_suffix
		log_info(f"{TAB * indent}Queueing {name}")

		queue_spec = {'job_name': name, 'cpus_per_task': cpus}
		spec: Dict[str, Any] = {'_queueadapter': queue_spec}
		if priority is not None:
			spec['_priority'] = priority

		firework = Firework(firetask, name=name, spec=spec, parents=parents)
		self.wf_fws.append(firework)
		return firework

	def add_links(self, parent: Firework, *children: Firework) -> None:
		"""Add parent -> child dependency links."""
		self.wf_links[parent].extend(children)

	def build_wcm(self, operons: str) -> Firework:
		"""Build a Whole Cell Model workflow (Parca, sims, analysis, and output
		file compression).
		"""
		self.operons = operons
		self.name_suffix = (
			constants.OPERON_SUFFIX if OPERONS == 'both' and operons == 'on' else '')

		log_info(f"\n--- Building a WCM workflow with {operons=} ---")
		self.make_output_directories()
		self.write_metadata()
		return self.build_wcm_firetasks()

	# noinspection PyUnusedLocal
	def make_output_directories(self) -> None:
		"""Make the output directories for the Firetasks and set self.* fields
		to some of those paths.
		"""
		# NOTE: Beyond setting the instance variables, making the dirs should be
		# superfluous now that these Firetasks make their output dirs to support
		# cloud workflows.
		log_info("Making the output directories.")

		self.INDIV_OUT_DIRECTORY = filepath.makedirs(
			OUT_DIRECTORY, f"{SUBMISSION_TIME}__{SIM_DESCRIPTION}{self.name_suffix}")
		self.KB_DIRECTORY = filepath.makedirs(self.INDIV_OUT_DIRECTORY, constants.KB_DIR)
		self.VARIANT_PLOT_DIRECTORY = os.path.join(self.INDIV_OUT_DIRECTORY, constants.PLOTOUT_DIR)

		for i in VARIANTS_TO_RUN:
			VARIANT_DIRECTORY = filepath.makedirs(self.INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
			VARIANT_SIM_DATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, constants.VKB_DIR)
			VARIANT_METADATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, constants.METADATA_DIR)
			VARIANT_COHORT_PLOT_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, constants.PLOTOUT_DIR)

			for j in range(SEED, SEED + N_INIT_SIMS):
				SEED_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "%06d" % j)
				SEED_PLOT_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, constants.PLOTOUT_DIR)

				for k in range(N_GENS):
					GEN_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, "generation_%06d" % k)

					for l in (range(2**k) if not SINGLE_DAUGHTERS else [0]):
						CELL_DIRECTORY = filepath.makedirs(GEN_DIRECTORY, "%06d" % l)
						CELL_SIM_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "simOut")
						CELL_PLOT_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, constants.PLOTOUT_DIR)
						CELL_SERIES_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "seriesOut")

	def write_metadata(self):
		"""Write the metadata files and set self.metadata."""
		self.metadata = {
			"git_hash": filepath.git_hash(),
			"git_branch": filepath.git_branch(),
			"description": os.environ.get("DESC", ""),
			"operons": self.operons,
			"time": SUBMISSION_TIME,
			"python": sys.version.splitlines()[0],
			"total_gens": N_GENS,
			"analysis_type": None,
			"variant": VARIANT,
			"total_variants": str(len(VARIANTS_TO_RUN)),
			"mass_distribution": MASS_DISTRIBUTION,
			"growth_rate_noise": GROWTH_RATE_NOISE,
			"d_period_division": D_PERIOD_DIVISION,
			"variable_elongation_transcription": VARIABLE_ELONGATION_TRANSCRIPTION,
			"variable_elongation_translation": VARIABLE_ELONGATION_TRANSLATION,
			"translation_supply": TRANSLATION_SUPPLY,
			"trna_charging": TRNA_CHARGING,
			"aa_supply_in_charging": AA_SUPPLY_IN_CHARGING,
			"ppgpp_regulation": PPGPP_REGULATION,
			"superhelical_density": SUPERHELICAL_DENSITY,
			"recycle_stalled_elongation": RECYCLE_STALLED_ELONGATION,
			"mechanistic_replisome": MECHANISTIC_REPLISOME,
			"mechanistic_translation_supply": MECHANISTIC_TRANSLATION_SUPPLY,
			"mechanistic_aa_transport": MECHANISTIC_AA_TRANSPORT,
			"trna_attenuation": TRNA_ATTENUATION,
			"adjust_timestep_for_charging": ADJUST_TIMESTEP_FOR_CHARGING,
			}

		METADATA_DIRECTORY = filepath.makedirs(self.INDIV_OUT_DIRECTORY, constants.METADATA_DIR)
		metadata_path = os.path.join(METADATA_DIRECTORY, constants.JSON_METADATA_FILE)
		filepath.write_json_file(metadata_path, self.metadata)

		git_diff = filepath.run_cmdline("git diff HEAD", trim=False)
		if git_diff:
			filepath.write_file(os.path.join(METADATA_DIRECTORY, "git_diff.txt"), git_diff)

	def build_wcm_firetasks(self):
		"""Build the WCM Firetasks and their dependency links.

		Call convert_to_fireworks_workflow() to convert the accumulated info to
		a Fireworks Workflow object.

		Returns a parent for comparison analysis (between this WCM and others)
		Firetasks **only if RUN_AGGREGATE_ANALYSIS**.

		**Priorities** are set to reduce possible worker starvation time
		(waiting on upstream tasks) by favoring tasks with more downstream tasks.
		Parca and variant: 12, sim 10-11, analysis 2-5, compression ≤ 0.
		"""
		INDIV_OUT_DIRECTORY = self.INDIV_OUT_DIRECTORY
		KB_DIRECTORY = self.KB_DIRECTORY
		VARIANT_PLOT_DIRECTORY = self.VARIANT_PLOT_DIRECTORY
		build_causality_network = BUILD_CAUSALITY_NETWORK

		# Initialize KB
		fw_init_raw_data = self.add_firework(
			InitRawDataTask(
				operons=self.operons,
				output=os.path.join(KB_DIRECTORY, constants.SERIALIZED_RAW_DATA)),
			"InitRawData",
			priority=12)

		# CalculateSimData
		cpusForParca = 8 if PARALLEL_PARCA else 1
		fw_calculate_sim_data = self.add_firework(
			FitSimDataTask(
				input_data=os.path.join(KB_DIRECTORY, constants.SERIALIZED_RAW_DATA),
				output_data=os.path.join(KB_DIRECTORY, constants.SERIALIZED_SIM_DATA_FILENAME),
				cached=CACHED_SIM_DATA,
				cached_data=os.path.join(CACHED_SIM_DATA_DIRECTORY,
										 constants.SERIALIZED_SIM_DATA_FILENAME),
				cpus=cpusForParca,
				debug=DEBUG_PARCA,
				disable_ribosome_capacity_fitting=DISABLE_RIBOSOME_CAPACITY_FITTING,
				disable_rnapoly_capacity_fitting=DISABLE_RNAPOLY_CAPACITY_FITTING,
				output_metrics_data=os.path.join(
					KB_DIRECTORY, constants.SERIALIZED_METRICS_DATA_FILENAME)),
			name="CalculateSimData",
			parents=fw_init_raw_data,
			cpus=cpusForParca,
			priority=12)

		# Raw KB compression
		fw_raw_data_compression = None
		if COMPRESS_OUTPUT:
			fw_raw_data_compression = self.add_firework(
				ScriptTask(
					script="bzip2 -v " + os.path.join(KB_DIRECTORY, constants.SERIALIZED_RAW_DATA)),
				name="ScriptTask_compression_raw_data",
				parents=fw_calculate_sim_data)

		# SimData compression
		fw_sim_data_1_compression = None
		if COMPRESS_OUTPUT:
			fw_sim_data_1_compression = self.add_firework(
				ScriptTask(
					script="bzip2 -v " + os.path.join(KB_DIRECTORY,
													  constants.SERIALIZED_SIM_DATA_FILENAME)),
				name="ScriptTask_compression_sim_data",
				parents=fw_calculate_sim_data)

		# Initiate raw validation data
		fw_raw_validation_data = self.add_firework(
			InitRawValidationDataTask(
				output=os.path.join(KB_DIRECTORY, constants.SERIALIZED_RAW_VALIDATION_DATA)),
			name="InitValidationDataRaw",
			priority=12)

		# Raw validation data compression
		fw_raw_validation_data_compression = None
		if COMPRESS_OUTPUT:
			fw_raw_validation_data_compression = self.add_firework(
				ScriptTask(
					script="bzip2 -v " + os.path.join(KB_DIRECTORY,
													  constants.SERIALIZED_RAW_VALIDATION_DATA)),
				name="ScriptTask_compression_validation_data_raw")


		# Initialize full validation data
		fw_validation_data = self.add_firework(
			InitValidationDataTask(
				validation_data_input=os.path.join(KB_DIRECTORY,
												   constants.SERIALIZED_RAW_VALIDATION_DATA),
				knowledge_base_raw=os.path.join(KB_DIRECTORY, constants.SERIALIZED_RAW_DATA),
				output_data=os.path.join(KB_DIRECTORY, constants.SERIALIZED_VALIDATION_DATA)),
			name="InitValidationData",
			parents=[fw_raw_validation_data, fw_init_raw_data],
			priority=12)

		# Full validation data compression
		fw_validation_data_compression = None
		if COMPRESS_OUTPUT:
			fw_validation_data_compression = self.add_firework(
				ScriptTask(
					script="bzip2 -v " + os.path.join(KB_DIRECTORY,
													  constants.SERIALIZED_VALIDATION_DATA)),
				name="ScriptTask_compression_validation_data")
			self.add_links(fw_validation_data,
						   fw_raw_validation_data_compression, fw_raw_data_compression)

		# Parca analysis
		fw_variant_analysis = None

		if RUN_AGGREGATE_ANALYSIS:
			fw_parca_analysis = self.add_firework(
				AnalysisParcaTask(
					input_directory=KB_DIRECTORY,
					input_sim_data=os.path.join(KB_DIRECTORY,
												constants.SERIALIZED_SIM_DATA_FILENAME),
					input_validation_data=os.path.join(KB_DIRECTORY,
													   constants.SERIALIZED_VALIDATION_DATA),
					output_plots_directory=os.path.join(INDIV_OUT_DIRECTORY,
														constants.KB_PLOT_OUTPUT_DIR),
					plot=PLOTS,
					cpus=analysis_cpus,
					metadata=self.metadata),
				name="AnalysisParcaTask",
				parents=[fw_calculate_sim_data, fw_validation_data],
				cpus=analysis_cpus,
				priority=5)
			if COMPRESS_OUTPUT:
				self.add_links(fw_parca_analysis,
							   fw_sim_data_1_compression, fw_validation_data_compression)

			# Variant analysis
			self.fw_variant_analysis = fw_variant_analysis = self.add_firework(
				AnalysisVariantTask(
					input_directory=INDIV_OUT_DIRECTORY,
					input_sim_data=os.path.join(KB_DIRECTORY,
												constants.SERIALIZED_SIM_DATA_FILENAME),
					input_validation_data=os.path.join(KB_DIRECTORY,
													   constants.SERIALIZED_VALIDATION_DATA),
					output_plots_directory=VARIANT_PLOT_DIRECTORY,
					plot=PLOTS,
					cpus=analysis_cpus,
					metadata=self.metadata),
				name="AnalysisVariantTask",
				cpus=analysis_cpus,
				priority=5)

		### Create variants and simulations
		fw_this_variant_sim_data_compression = None
		fw_this_variant_this_gen_this_sim_compression = None

		for i in VARIANTS_TO_RUN:
			VARIANT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
			VARIANT_SIM_DATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, constants.VKB_DIR)
			VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, constants.METADATA_DIR)
			md_cohort = dict(self.metadata, variant_function = VARIANT, variant_index = i)

			# Variant simData creation task
			# Note: This task doesn't depend on fw_validation_data but such a
			# link is lighter weight for Fireworks than making every analysis
			# task directly depend on it.
			fw_this_variant_sim_data = self.add_firework(
				VariantSimDataTask(
					variant_function=VARIANT,
					variant_index=i,
					input_sim_data=os.path.join(KB_DIRECTORY,
												constants.SERIALIZED_SIM_DATA_FILENAME),
					output_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
												 constants.SERIALIZED_SIM_DATA_MODIFIED),
					variant_metadata_directory=VARIANT_METADATA_DIRECTORY),
				name=f"VariantSimDataTask__{VARIANT}_{i:06d}",
				parents=[fw_calculate_sim_data, fw_validation_data],
				priority=12)

			if COMPRESS_OUTPUT:
				self.add_links(fw_this_variant_sim_data, fw_sim_data_1_compression)

				# Variant simData compression
				fw_this_variant_sim_data_compression = self.add_firework(
					ScriptTask(
						script="bzip2 -v " + os.path.join(VARIANT_SIM_DATA_DIRECTORY,
														  constants.SERIALIZED_SIM_DATA_MODIFIED)),
					name="ScriptTask_compression_variant_KB")

			# Cohort analysis
			COHORT_PLOT_DIRECTORY = os.path.join(VARIANT_DIRECTORY, constants.PLOTOUT_DIR)

			fw_this_variant_cohort_analysis = None

			if RUN_AGGREGATE_ANALYSIS:
				fw_this_variant_cohort_analysis = self.add_firework(
					AnalysisCohortTask(
						input_variant_directory=VARIANT_DIRECTORY,
						input_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
													constants.SERIALIZED_SIM_DATA_MODIFIED),
						input_validation_data=os.path.join(KB_DIRECTORY,
														   constants.SERIALIZED_VALIDATION_DATA),
						output_plots_directory=COHORT_PLOT_DIRECTORY,
						plot=PLOTS,
						cpus=analysis_cpus,
						metadata=md_cohort),
					name=f"AnalysisCohortTask__Var_{i:02d}",
					cpus=analysis_cpus,
					priority=4)

			fw_this_variant_this_seed_multigen_analysis = None

			for j in range(SEED, SEED + N_INIT_SIMS):
				SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
				SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, constants.PLOTOUT_DIR)
				md_multigen = dict(md_cohort, seed = j)

				if RUN_AGGREGATE_ANALYSIS:
					fw_this_variant_this_seed_multigen_analysis = self.add_firework(
						AnalysisMultiGenTask(
							input_seed_directory=SEED_DIRECTORY,
							input_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
														constants.SERIALIZED_SIM_DATA_MODIFIED),
							input_validation_data=os.path.join(KB_DIRECTORY,
															   constants.SERIALIZED_VALIDATION_DATA),
							output_plots_directory=SEED_PLOT_DIRECTORY,
							plot=PLOTS,
							cpus=analysis_cpus,
							metadata=md_multigen),
						name=f"AnalysisMultiGenTask__Var_{i:02d}__Seed_{j:06d}",
						cpus=analysis_cpus,
						priority=3,
						indent=1)

					if COMPRESS_OUTPUT:
						self.add_links(fw_this_variant_this_seed_multigen_analysis,
									   fw_this_variant_sim_data_compression)

				sims_this_seed = collections.defaultdict(list)

				for k in range(N_GENS):
					GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)
					md_single = dict(md_multigen, gen = k)

					for l in (range(2**k) if not SINGLE_DAUGHTERS else [0]):
						CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
						CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
						CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, constants.PLOTOUT_DIR)
						CELL_SERIES_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "seriesOut")

						# Simulation task
						sim_fw_name = f"SimulationTask__Var_{i:02d}__Seed_{j:d}__Gen_{k:d}__Cell_{l:d}"
						sim_task_args = dict(
							input_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
														constants.SERIALIZED_SIM_DATA_MODIFIED),
							output_directory=CELL_SIM_OUT_DIRECTORY,
							timeline=TIMELINE,
							length_sec=WC_LENGTHSEC,
							timestep_safety_frac=TIMESTEP_SAFETY_FRAC,
							timestep_max=TIMESTEP_MAX,
							timestep_update_freq=TIMESTEP_UPDATE_FREQ,
							adjust_timestep_for_charging=ADJUST_TIMESTEP_FOR_CHARGING,
							log_to_disk_every=LOG_TO_DISK_EVERY,
							jit=JIT,
							mass_distribution=MASS_DISTRIBUTION,
							growth_rate_noise=GROWTH_RATE_NOISE,
							d_period_division=D_PERIOD_DIVISION,
							variable_elongation_transcription=VARIABLE_ELONGATION_TRANSCRIPTION,
							variable_elongation_translation=VARIABLE_ELONGATION_TRANSLATION,
							translation_supply=TRANSLATION_SUPPLY,
							trna_charging=TRNA_CHARGING,
							aa_supply_in_charging=AA_SUPPLY_IN_CHARGING,
							ppgpp_regulation=PPGPP_REGULATION,
							superhelical_density=SUPERHELICAL_DENSITY,
							recycle_stalled_elongation=RECYCLE_STALLED_ELONGATION,
							mechanistic_replisome=MECHANISTIC_REPLISOME,
							mechanistic_translation_supply=MECHANISTIC_TRANSLATION_SUPPLY,
							mechanistic_aa_transport=MECHANISTIC_AA_TRANSPORT,
							trna_attenuation=TRNA_ATTENUATION,
							raise_on_time_limit=RAISE_ON_TIME_LIMIT)

						if k == 0:
							fw_this_variant_this_gen_this_sim = self.add_firework(
								SimulationTask(
									seed=j,
									**sim_task_args),
								name=sim_fw_name,
								cpus=1,
								priority=10,
								indent=2)
						elif k > 0:
							PARENT_GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % (k - 1))
							PARENT_CELL_DIRECTORY = os.path.join(PARENT_GEN_DIRECTORY, "%06d" % (l // 2))
							PARENT_CELL_SIM_OUT_DIRECTORY = os.path.join(PARENT_CELL_DIRECTORY, "simOut")
							DAUGHTER_STATE_PATH = os.path.join(PARENT_CELL_SIM_OUT_DIRECTORY,
								constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))

							fw_this_variant_this_gen_this_sim = self.add_firework(
								SimulationDaughterTask(
									inherited_state_path=DAUGHTER_STATE_PATH,
									seed=(j + 1) * ((2 ** k - 1) + l),
									**sim_task_args),
								name=sim_fw_name,
								cpus=1,
								priority=11,
								indent=2)
						else:
							raise ValueError("k ({}) < 0".format(k))

						# Only add the last generation as dependencies for multiple sim analysis tasks
						if RUN_AGGREGATE_ANALYSIS and k == N_GENS - 1:
							self.add_links(fw_this_variant_this_gen_this_sim,
										   fw_this_variant_this_seed_multigen_analysis,
										   fw_this_variant_cohort_analysis,
										   fw_variant_analysis)

						sims_this_seed[k].append(fw_this_variant_this_gen_this_sim)

						if k == 0:
							self.add_links(fw_this_variant_sim_data, fw_this_variant_this_gen_this_sim)

						elif k > 0:
							fw_parent_sim = sims_this_seed[k - 1][l // 2]
							self.add_links(fw_parent_sim, fw_this_variant_this_gen_this_sim)

						if COMPRESS_OUTPUT:
							# Output compression job
							fw_this_variant_this_gen_this_sim_compression = self.add_firework(
								ScriptTask(
									script='for dir in %s; do echo "Compressing $dir"; find "$dir" -type f | xargs bzip2; done' % os.path.join(
										CELL_SIM_OUT_DIRECTORY, "*")),
								name=f"ScriptTask_compression_simulation__Seed_{j:d}__Gen_{k:d}__Cell_{l:d}",
								priority=0)

						if RUN_AGGREGATE_ANALYSIS:
							# AnalysisSingle task
							fw_this_variant_this_gen_this_sim_analysis = self.add_firework(
								AnalysisSingleTask(
									input_results_directory=CELL_SIM_OUT_DIRECTORY,
									input_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
																constants.SERIALIZED_SIM_DATA_MODIFIED),
									input_validation_data=os.path.join(KB_DIRECTORY,
																	   constants.SERIALIZED_VALIDATION_DATA),
									output_plots_directory=CELL_PLOT_OUT_DIRECTORY,
									plot=PLOTS,
									cpus=analysis_cpus,
									metadata=md_single),
								name=f"AnalysisSingleTask__Var_{i:d}__Seed_{j:d}__Gen_{k:d}__Cell_{l:d}",
								parents=fw_this_variant_this_gen_this_sim,
								cpus=analysis_cpus,
								priority=2,
								indent=3)

							if COMPRESS_OUTPUT:
								# Don't compress any outputs or validation data until all analysis scripts (single gen, multigen, and cohort) have finished running
								data_fws = [
									fw_this_variant_this_gen_this_sim_analysis,
									fw_this_variant_this_seed_multigen_analysis,
									fw_this_variant_cohort_analysis,
									fw_variant_analysis]
								for data in data_fws:
									self.add_links(data,
										fw_this_variant_sim_data_compression,
										fw_validation_data_compression,
										fw_this_variant_this_gen_this_sim_compression)


						if build_causality_network:
							# BuildCausalityNetwork task
							build_causality_network = False  # once per WCM is enough
							fw_this_variant_this_gen_this_sim_causality_network = self.add_firework(
								BuildCausalityNetworkTask(
									input_results_directory=CELL_SIM_OUT_DIRECTORY,
									input_sim_data=os.path.join(VARIANT_SIM_DATA_DIRECTORY,
																constants.SERIALIZED_SIM_DATA_MODIFIED),
									output_dynamics_directory=CELL_SERIES_OUT_DIRECTORY,
									metadata=md_single,
									),
								name=f"BuildCausalityNetworkTask__Var_{i:d}__Seed_{j:d}__Gen_{k:d}__Cell_{l:d}",
								parents=fw_this_variant_this_gen_this_sim,
								cpus=analysis_cpus,
								priority=2,
								indent=3)

							if COMPRESS_OUTPUT:
								# Don't compress any outputs or sim_data until
								# causality network scripts have finished running
								self.add_links(fw_this_variant_this_gen_this_sim_causality_network,
											   fw_this_variant_sim_data_compression,
											   fw_this_variant_this_gen_this_sim_compression)

	def add_comparison_analysis(self, reference_sim_dir: str,
			reference_variant_analysis: Firework) -> None:
		"""Add an AnalysisComparisonTask that compares the WCM workflow results
		against a reference WCM workflow's results.
		"""
		if RUN_AGGREGATE_ANALYSIS:
			plot_out_dir = os.path.join(
				self.INDIV_OUT_DIRECTORY, constants.COMPARISON_PLOTOUT_DIR)
			self.add_firework(
				AnalysisComparisonTask(
					reference_sim_dir=reference_sim_dir,
					input_sim_dir=self.INDIV_OUT_DIRECTORY,
					output_plots_directory=plot_out_dir,
					metadata=self.metadata,
					plot=PLOTS,
					cpus=analysis_cpus),
				name="AnalysisComparisonTask",
				parents=[reference_variant_analysis, self.fw_variant_analysis],
				cpus=analysis_cpus,
				priority=2)

	def convert_to_fireworks_workflow(self) -> Workflow:
		"""Build a Fireworks Workflow object from the Firetask and dependency
		info accumulated by build_wcm().
		"""
		workflow = Workflow(self.wf_fws, links_dict=self.wf_links)
		return workflow


def upload_workflow(workflow: Workflow):
	"""Upload a Fireworks Workflow to the LaunchPad."""
	with open(LAUNCHPAD_FILE) as f:
		lpad = LaunchPad(**yaml.safe_load(f))
	lpad.add_wf(workflow)


def main():
	builder = WorkflowBuilder()

	if OPERONS == 'both':
		builder.build_wcm(operons='off')
		sim_dir1 = builder.INDIV_OUT_DIRECTORY
		variant_analysis1 = builder.fw_variant_analysis

		builder.build_wcm(operons='on')

		builder.add_comparison_analysis(sim_dir1, variant_analysis1)
	else:
		builder.build_wcm(operons=OPERONS)

	wf = builder.convert_to_fireworks_workflow()
	# (Could optionally dump wf as a YAML file here.)
	upload_workflow(wf)


if __name__ == '__main__':
	main()
