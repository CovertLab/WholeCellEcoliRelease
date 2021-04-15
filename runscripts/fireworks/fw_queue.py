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
	RUN_AGGREGATE_ANALYSIS (int, "1"): if nonzero, all analyses are run on
		simulation output
	PLOTS (str, "CORE"): Which analyses to run (if RUN_AGGREGATE_ANALYSIS is
		true). This should name one or more tags. For more than one tag,
		separate them with whitespace and remember shell quoting. ACTIVE
		includes all active plots. CORE includes just the plots recommended for
		everyday development. You can also name specific analysis files but any
		analysis categories that don't have such a filename will print error
		messages.
	DISABLE_RIBOSOME_CAPACITY_FITTING (int, "0"): if nonzero, ribosome
		expression is not fit to protein synthesis demands
	DISABLE_RNAPOLY_CAPACITY_FITTING (int, "0"): if nonzero, RNA polymerase
		expression is not fit to RNA synthesis demands
	WC_ANALYZE_FAST (anything, --): if set, run each analysis plot in a separate
		process
	BUILD_CAUSALITY_NETWORK (int, "0"): if nonzero, causality network files are
		generated from simulation output
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
	TRANSLATION_SUPPLY (int, "1"): if nonzero, the ribosome elongation rate is
		limited by the condition specific rate of amino acid supply; otherwise
		the elongation rate is set by condition
	TRNA_CHARGING (int, "1"): if nonzero, tRNA charging reactions are modeled
		and the ribosome elongation rate is set by the amount of charged tRNA
		present.  This option will override TRANSLATION_SUPPLY in the simulation.
	PPGPP_REGULATION (int, "0"): if nonzero, ppGpp concentration is determined
		with kinetic equations
	SUPERHELICAL_DENSITY (int, "0"): if nonzero, dynamically compute
		superhelical densities of each DNA fragment
	RECYCLE_STALLED_ELONGATION (int "0"): if nonzero, recycle RNAP and fragment
		bases when transcription elongation is stalled in ntp-limiting conditions
	MECHANISTIC_REPLISOME (int, "1"): if nonzero, replisome initiation is
		mechanistic (requires appropriate number of subunits to initiate)
	MECHANISTIC_AA_SUPPLY (int, "0"): if nonzero, amino acid supply is
		mechanistic (depends on concentrations of enzymes and amino acids)
	TRNA_ATTENUATION (int, "0"): if nonzero, transcriptional attenuation by
		charged tRNA is enabled

Additional variables:
	LAUNCHPAD_FILE (str, "my_launchpad.yaml"): set launchpad config file location
	VERBOSE_QUEUE (int, "1"): if nonzero, gives more detailed messages during
		fireworks set up

Environment variables that matter when running the workflow:
	DEBUG_GC (int, "0"): if nonzero, enable leak detection in the analysis plots


-------------------------------------------------
Current dependency network of existing FireTasks (each task has a list of its
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

FitSimData (fw_calculate_sim_data)
* CompressRawData if COMPRESS_OUTPUT
* VariantSimData * VARIANTS_TO_RUN
* AnalysisParca if RUN_AGGREGATE_ANALYSIS

VariantSimData (fw_this_variant_sim_data)
* CompressFitSimData if COMPRESS_OUTPUT
* Simulation/SimulationDaughter * N_INIT_SIMS

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
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

AnalysisCohort (fw_this_variant_cohort_analysis)
* CompressVariantSimData if COMPRESS_OUTPUT
* CompressValidationData if COMPRESS_OUTPUT
* CompressSimulationOutput if COMPRESS_OUTPUT

AnalysisMultiGen (fw_this_variant_this_seed_this_analysis)
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

from __future__ import absolute_import, division, print_function

import collections
import os
import sys

import yaml
from fireworks import Firework, LaunchPad, Workflow, ScriptTask

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
from wholecell.fireworks.firetasks import BuildCausalityNetworkTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from wholecell.utils import constants
from wholecell.utils import filepath
from six.moves import range


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

print('Parsed environmental variables:')

#### Initial setup ###

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
LOG_TO_DISK_EVERY = int(get_environment("LOG_TO_DISK_EVERY", DEFAULT_SIMULATION_KWARGS["logToDiskEvery"]))
JIT = bool(int(get_environment("JIT", DEFAULT_SIMULATION_KWARGS["jit"])))
MASS_DISTRIBUTION = bool(int(get_environment("MASS_DISTRIBUTION", DEFAULT_SIMULATION_KWARGS["massDistribution"])))
GROWTH_RATE_NOISE = bool(int(get_environment("GROWTH_RATE_NOISE", DEFAULT_SIMULATION_KWARGS["growthRateNoise"])))
D_PERIOD_DIVISION = bool(int(get_environment("D_PERIOD_DIVISION", DEFAULT_SIMULATION_KWARGS["dPeriodDivision"])))
TRANSLATION_SUPPLY = bool(int(get_environment("TRANSLATION_SUPPLY", DEFAULT_SIMULATION_KWARGS["translationSupply"])))
TRNA_CHARGING = bool(int(get_environment("TRNA_CHARGING", DEFAULT_SIMULATION_KWARGS["trna_charging"])))
PPGPP_REGULATION = bool(int(get_environment("PPGPP_REGULATION", DEFAULT_SIMULATION_KWARGS["ppgpp_regulation"])))
SUPERHELICAL_DENSITY = bool(int(get_environment("SUPERHELICAL_DENSITY", DEFAULT_SIMULATION_KWARGS["superhelical_density"])))
RECYCLE_STALLED_ELONGATION = bool(int(get_environment("RECYCLE_STALLED_ELONGATION", DEFAULT_SIMULATION_KWARGS["recycle_stalled_elongation"])))
MECHANISTIC_REPLISOME = bool(int(get_environment("MECHANISTIC_REPLISOME", DEFAULT_SIMULATION_KWARGS["mechanistic_replisome"])))
MECHANISTIC_AA_SUPPLY = bool(int(get_environment("MECHANISTIC_AA_SUPPLY", DEFAULT_SIMULATION_KWARGS["mechanistic_aa_supply"])))
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
PLOTS = get_environment("PLOTS", "CORE").split()
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


### Set path variables and create directories

OUT_DIRECTORY = filepath.makedirs(filepath.ROOT_PATH, "out")
CACHED_SIM_DATA_DIRECTORY = os.path.join(filepath.ROOT_PATH, "cached")

# To run each analysis plot in a separate process, ask the analysis Firetasks
# for several CPUs (it will clip to the number available) and allocate multiple
# CPUs from SLURM via the Fireworks Queue Adaptor. Otherwise, let the qadapter
# YAML file request the number of CPUs so we can tune it to ask for extra CPUs
# in order to get proportionally more RAM, e.g. after running many generations.
if "WC_ANALYZE_FAST" in os.environ:
	analysis_cpus = 8
	analysis_q_cpus = {"cpus_per_task": analysis_cpus}
else:
	analysis_cpus = 1
	analysis_q_cpus = {}

SUBMISSION_TIME = filepath.timestamp()
INDIV_OUT_DIRECTORY = filepath.makedirs(OUT_DIRECTORY, SUBMISSION_TIME + "__" + SIM_DESCRIPTION)
KB_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, "kb")
METADATA_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, "metadata")


log_info("Building filestructure.")

for i in VARIANTS_TO_RUN:
	VARIANT_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "metadata")
	VARIANT_COHORT_PLOT_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "plotOut")

	for j in range(SEED, SEED + N_INIT_SIMS):
		SEED_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, "plotOut")

		for k in range(N_GENS):
			GEN_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, "generation_%06d" % k)

			for l in (range(2**k) if not SINGLE_DAUGHTERS else [0]):
				CELL_DIRECTORY = filepath.makedirs(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "plotOut")
				CELL_SERIES_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "seriesOut")

### Write metadata
metadata = {
	"git_hash": filepath.git_hash(),
	"git_branch": filepath.git_branch(),
	"description": os.environ.get("DESC", ""),
	"time": SUBMISSION_TIME,
	"python": sys.version.splitlines()[0],
	"total_gens": N_GENS,
	"analysis_type": None,
	"variant": VARIANT,
	"total_variants": str(len(VARIANTS_TO_RUN)),
	"mass_distribution": MASS_DISTRIBUTION,
	"growth_rate_noise": GROWTH_RATE_NOISE,
	"d_period_division": D_PERIOD_DIVISION,
	"translation_supply": TRANSLATION_SUPPLY,
	"trna_charging": TRNA_CHARGING,
	"ppgpp_regulation": PPGPP_REGULATION,
	"superhelical_density": SUPERHELICAL_DENSITY,
	"recycle_stalled_elongation": RECYCLE_STALLED_ELONGATION,
	"mechanistic_replisome": MECHANISTIC_REPLISOME,
	"mechanistic_aa_supply": MECHANISTIC_AA_SUPPLY,
	"trna_attenuation": TRNA_ATTENUATION,
	}

metadata_path = os.path.join(METADATA_DIRECTORY, constants.JSON_METADATA_FILE)
filepath.write_json_file(metadata_path, metadata)

git_diff = filepath.run_cmdline("git diff HEAD", trim=False)
if git_diff:
	filepath.write_file(os.path.join(METADATA_DIRECTORY, "git_diff.txt"), git_diff)

#### Create workflow

# Create launchpad
with open(LAUNCHPAD_FILE) as f:
	lpad = LaunchPad(**yaml.safe_load(f))

# Store list of FireWorks
wf_fws = []

# Store links defining parent/child dependency relationships of FireWorks
wf_links = collections.defaultdict(list)


### Initialize KB

filename_raw_data = constants.SERIALIZED_RAW_DATA
filename_sim_data = constants.SERIALIZED_SIM_DATA_FILENAME
filename_sim_data_modified = constants.SERIALIZED_SIM_DATA_MODIFIED

fw_name = "InitRawData"

log_info("Queueing {}".format(fw_name))

fw_init_raw_data = Firework(
	InitRawDataTask(
		output = os.path.join(KB_DIRECTORY, filename_raw_data)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_init_raw_data)

### Fit (Level 1)

fw_name = "CalculateSimData"

log_info("Queueing {}".format(fw_name))

cpusForParca = 8 if PARALLEL_PARCA else 1
fw_calculate_sim_data = Firework(
	FitSimDataTask(
		input_data = os.path.join(KB_DIRECTORY, filename_raw_data),
		output_data = os.path.join(KB_DIRECTORY, filename_sim_data),
		cached = CACHED_SIM_DATA,
		cached_data = os.path.join(CACHED_SIM_DATA_DIRECTORY, filename_sim_data),
		cpus = cpusForParca,
		debug = DEBUG_PARCA,
		disable_ribosome_capacity_fitting = DISABLE_RIBOSOME_CAPACITY_FITTING,
		disable_rnapoly_capacity_fitting = DISABLE_RNAPOLY_CAPACITY_FITTING,
		output_metrics_data = os.path.join(
			KB_DIRECTORY, constants.SERIALIZED_METRICS_DATA_FILENAME),
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": cpusForParca}, "_priority":1}
	)

wf_fws.append(fw_calculate_sim_data)
wf_links[fw_init_raw_data].append(fw_calculate_sim_data)

# Unfit KB compression
fw_raw_data_compression = None
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_raw_data"

	log_info("Queueing {}".format(fw_name))

	fw_raw_data_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_raw_data)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_raw_data_compression)
	wf_links[fw_calculate_sim_data].append(fw_raw_data_compression)

# Fit Level 1 KB compression

fw_sim_data_1_compression = None
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_sim_data"

	log_info("Queueing {}".format(fw_name))

	fw_sim_data_1_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_sim_data)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_sim_data_1_compression)


### Initialize validation data

# Initiate raw validation data
filename_raw_validation_data = constants.SERIALIZED_RAW_VALIDATION_DATA

fw_name = "InitValidationDataRaw"

log_info("Queueing {}".format(fw_name))

fw_raw_validation_data = Firework(
	InitRawValidationDataTask(
		output = os.path.join(KB_DIRECTORY, filename_raw_validation_data)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_raw_validation_data)

# Raw validation data compression
fw_raw_validation_data_compression = None
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_validation_data_raw"
	fw_raw_validation_data_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_raw_validation_data)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_raw_validation_data_compression)


# Initialize full validation data
filename_validation_data = constants.SERIALIZED_VALIDATION_DATA

fw_name = "InitValidationData"

log_info("Queueing {}".format(fw_name))

fw_validation_data = Firework(
	InitValidationDataTask(
		validation_data_input = os.path.join(KB_DIRECTORY, filename_raw_validation_data),
		knowledge_base_raw = os.path.join(KB_DIRECTORY, filename_raw_data),
		output_data = os.path.join(KB_DIRECTORY, filename_validation_data),
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_validation_data)
wf_links[fw_raw_validation_data].append(fw_validation_data)
wf_links[fw_init_raw_data].append(fw_validation_data)
fw_validation_data_compression = None

# Full validation data compression
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_validation_data"

	log_info("Queueing {}".format(fw_name))

	fw_validation_data_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_validation_data)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_validation_data_compression)

	wf_links[fw_validation_data].append(fw_raw_validation_data_compression)
	wf_links[fw_validation_data].append(fw_raw_data_compression)

# Variant analysis
VARIANT_PLOT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, "plotOut")
fw_variant_analysis = None
fw_parca_analysis = None

if RUN_AGGREGATE_ANALYSIS:
	fw_name = "AnalysisParcaTask"
	fw_parca_analysis = Firework(
		AnalysisParcaTask(
			input_directory = KB_DIRECTORY,
			input_sim_data = os.path.join(KB_DIRECTORY, filename_sim_data),
			input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
			output_plots_directory = os.path.join(INDIV_OUT_DIRECTORY, constants.KB_PLOT_OUTPUT_DIR),
			plot = PLOTS,
			cpus = analysis_cpus,
			metadata = metadata,
			),
		name = fw_name,
		spec = {"_queueadapter": dict(analysis_q_cpus, job_name=fw_name), "_priority":5}
		)
	wf_fws.append(fw_parca_analysis)
	wf_links[fw_calculate_sim_data].append(fw_parca_analysis)
	wf_links[fw_validation_data].append(fw_parca_analysis)
	if COMPRESS_OUTPUT:
		wf_links[fw_parca_analysis].append(fw_sim_data_1_compression)
		wf_links[fw_parca_analysis].append(fw_validation_data_compression)

	fw_name = "AnalysisVariantTask"
	fw_variant_analysis = Firework(
		AnalysisVariantTask(
			input_directory = os.path.join(INDIV_OUT_DIRECTORY),
			input_sim_data = os.path.join(KB_DIRECTORY, filename_sim_data),
			input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
			output_plots_directory = VARIANT_PLOT_DIRECTORY,
			plot = PLOTS,
			cpus = analysis_cpus,
			metadata = metadata,
			),
		name = fw_name,
		spec = {"_queueadapter": dict(analysis_q_cpus, job_name=fw_name), "_priority":5}
		)
	wf_fws.append(fw_variant_analysis)

### Create variants and simulations
fw_this_variant_sim_data_compression = None
fw_this_variant_this_gen_this_sim_compression = None

for i in VARIANTS_TO_RUN:
	log_info("Queueing Variant {} {}".format(VARIANT, i))
	VARIANT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")
	md_cohort = dict(metadata, variant_function = VARIANT, variant_index = i)

	# Variant simData creation task
	fw_name = "VariantSimDataTask_%06d" % (i,)
	fw_this_variant_sim_data = Firework(
		VariantSimDataTask(
			variant_function = VARIANT,
			variant_index = i,
			input_sim_data = os.path.join(KB_DIRECTORY, filename_sim_data),
			output_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
			variant_metadata_directory = VARIANT_METADATA_DIRECTORY,
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
		)

	wf_fws.append(fw_this_variant_sim_data)
	wf_links[fw_calculate_sim_data].append(fw_this_variant_sim_data)

	if COMPRESS_OUTPUT:
		wf_links[fw_this_variant_sim_data].append(fw_sim_data_1_compression)

		# Variant simData compression
		fw_name = "ScriptTask_compression_variant_KB"
		fw_this_variant_sim_data_compression = Firework(
			ScriptTask(
				script = "bzip2 -v " + os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified)
				),
			name = fw_name,
			spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
			)

		wf_fws.append(fw_this_variant_sim_data_compression)

	# Cohort analysis
	COHORT_PLOT_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "plotOut")

	fw_this_variant_cohort_analysis = None

	if RUN_AGGREGATE_ANALYSIS:
		fw_name = "AnalysisCohortTask__Var_%02d" % (i,)
		fw_this_variant_cohort_analysis = Firework(
			AnalysisCohortTask(
				input_variant_directory = VARIANT_DIRECTORY,
				input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
				input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
				output_plots_directory = COHORT_PLOT_DIRECTORY,
				plot = PLOTS,
				cpus = analysis_cpus,
				metadata = md_cohort,
				),
			name = fw_name,
			spec = {"_queueadapter": dict(analysis_q_cpus, job_name=fw_name), "_priority":4}
			)
		wf_fws.append(fw_this_variant_cohort_analysis)

	fw_this_variant_this_seed_this_analysis = None

	for j in range(SEED, SEED + N_INIT_SIMS):
		log_info("\tQueueing Seed {}".format(j))
		SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, "plotOut")
		md_multigen = dict(md_cohort, seed = j)

		if RUN_AGGREGATE_ANALYSIS:
			fw_name = "AnalysisMultiGenTask__Var_%02d__Seed_%06d" % (i, j)
			fw_this_variant_this_seed_this_analysis = Firework(
				AnalysisMultiGenTask(
					input_seed_directory = SEED_DIRECTORY,
					input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
					input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
					output_plots_directory = SEED_PLOT_DIRECTORY,
					plot = PLOTS,
					cpus = analysis_cpus,
					metadata = md_multigen,
					),
				name = fw_name,
				spec = {"_queueadapter": dict(analysis_q_cpus, job_name=fw_name), "_priority":3}
				)
			wf_fws.append(fw_this_variant_this_seed_this_analysis)

			if COMPRESS_OUTPUT:
				wf_links[fw_this_variant_this_seed_this_analysis].append(fw_this_variant_sim_data_compression)

		sims_this_seed = collections.defaultdict(list)

		for k in range(N_GENS):
			log_info("\t\tQueueing Gen %02d." % (k,))
			GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)
			md_single = dict(md_multigen, gen = k)

			for l in (range(2**k) if not SINGLE_DAUGHTERS else [0]):

				log_info("\t\t\tQueueing Cell {}".format(l))
				CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "plotOut")
				CELL_SERIES_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "seriesOut")

				# TODO: Add conditional logic here for mother vs daughter cells
				# Simulation task
				fw_name = "SimulationTask__Var_%02d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)

				if k == 0:
					fw_this_variant_this_gen_this_sim = Firework(
						SimulationTask(
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							seed = j,
							timeline = TIMELINE,
							length_sec = WC_LENGTHSEC,
							timestep_safety_frac = TIMESTEP_SAFETY_FRAC,
							timestep_max = TIMESTEP_MAX,
							timestep_update_freq = TIMESTEP_UPDATE_FREQ,
							log_to_disk_every = LOG_TO_DISK_EVERY,
							jit=JIT,
							mass_distribution = MASS_DISTRIBUTION,
							growth_rate_noise = GROWTH_RATE_NOISE,
							d_period_division = D_PERIOD_DIVISION,
							translation_supply = TRANSLATION_SUPPLY,
							trna_charging = TRNA_CHARGING,
							ppgpp_regulation = PPGPP_REGULATION,
							superhelical_density = SUPERHELICAL_DENSITY,
							recycle_stalled_elongation = RECYCLE_STALLED_ELONGATION,
							mechanistic_replisome = MECHANISTIC_REPLISOME,
							mechanistic_aa_supply = MECHANISTIC_AA_SUPPLY,
							trna_attenuation = TRNA_ATTENUATION,
							raise_on_time_limit = RAISE_ON_TIME_LIMIT,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": 1}, "_priority":10}
						)
				elif k > 0:
					PARENT_GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % (k - 1))
					PARENT_CELL_DIRECTORY = os.path.join(PARENT_GEN_DIRECTORY, "%06d" % (l // 2))
					PARENT_CELL_SIM_OUT_DIRECTORY = os.path.join(PARENT_CELL_DIRECTORY, "simOut")
					DAUGHTER_STATE_PATH = os.path.join(PARENT_CELL_SIM_OUT_DIRECTORY,
						constants.SERIALIZED_INHERITED_STATE % (l % 2 + 1))

					fw_this_variant_this_gen_this_sim = Firework(
						SimulationDaughterTask(
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							inherited_state_path = DAUGHTER_STATE_PATH,
							seed = (j + 1) * ((2**k - 1) + l),
							timeline = TIMELINE,
							length_sec = WC_LENGTHSEC,
							timestep_safety_frac = TIMESTEP_SAFETY_FRAC,
							timestep_max = TIMESTEP_MAX,
							timestep_update_freq = TIMESTEP_UPDATE_FREQ,
							log_to_disk_every = LOG_TO_DISK_EVERY,
							jit=JIT,
							mass_distribution = MASS_DISTRIBUTION,
							growth_rate_noise = GROWTH_RATE_NOISE,
							d_period_division = D_PERIOD_DIVISION,
							translation_supply = TRANSLATION_SUPPLY,
							trna_charging = TRNA_CHARGING,
							ppgpp_regulation = PPGPP_REGULATION,
							superhelical_density = SUPERHELICAL_DENSITY,
							recycle_stalled_elongation = RECYCLE_STALLED_ELONGATION,
							mechanistic_replisome = MECHANISTIC_REPLISOME,
							mechanistic_aa_supply = MECHANISTIC_AA_SUPPLY,
							trna_attenuation = TRNA_ATTENUATION,
							raise_on_time_limit = RAISE_ON_TIME_LIMIT,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": 1}, "_priority":11}
						)
				else:
					raise ValueError("k ({}) < 0".format(k))

				wf_fws.append(fw_this_variant_this_gen_this_sim)
				# Only add the last generation as dependencies for multiple sim analysis tasks
				if RUN_AGGREGATE_ANALYSIS and k == N_GENS - 1:
					wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_this_seed_this_analysis)
					wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_cohort_analysis)
					wf_links[fw_this_variant_this_gen_this_sim].append(fw_variant_analysis)

				sims_this_seed[k].append(fw_this_variant_this_gen_this_sim)

				if k == 0:
					wf_links[fw_this_variant_sim_data].append(fw_this_variant_this_gen_this_sim)

				elif k > 0:
					fw_parent_sim = sims_this_seed[k - 1][l // 2]
					wf_links[fw_parent_sim].append(fw_this_variant_this_gen_this_sim)

				if COMPRESS_OUTPUT:
					# Output compression job
					fw_name = "ScriptTask_compression_simulation__Seed_%d__Gen_%d__Cell_%d" % (
					j, k, l)
					fw_this_variant_this_gen_this_sim_compression = Firework(
						ScriptTask(
							script='for dir in %s; do echo "Compressing $dir"; find "$dir" -type f | xargs bzip2; done' % os.path.join(
								CELL_SIM_OUT_DIRECTORY, "*")
							),
						name=fw_name,
						spec={"_queueadapter": {"job_name": fw_name},
							"_priority": 0}
						)

					wf_fws.append(
						fw_this_variant_this_gen_this_sim_compression)

				if RUN_AGGREGATE_ANALYSIS:
					# AnalysisSingle task
					fw_name = "AnalysisSingleTask__Var_%d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)
					fw_this_variant_this_gen_this_sim_analysis = Firework(
						AnalysisSingleTask(
							input_results_directory = CELL_SIM_OUT_DIRECTORY,
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
							output_plots_directory = CELL_PLOT_OUT_DIRECTORY,
							plot = PLOTS,
							cpus = analysis_cpus,
							metadata = md_single,
							),
						name = fw_name,
						spec = {"_queueadapter": dict(analysis_q_cpus, job_name=fw_name), "_priority":2}
						)

					wf_fws.append(fw_this_variant_this_gen_this_sim_analysis)

					wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_this_gen_this_sim_analysis)

					if COMPRESS_OUTPUT:
						# Don't compress any outputs or validation data until all analysis scripts (single gen, multigen, and cohort) have finished running
						compression_fws = [
							fw_this_variant_sim_data_compression,
							fw_validation_data_compression,
							fw_this_variant_this_gen_this_sim_compression
						]
						data_fws = [
							fw_this_variant_this_gen_this_sim_analysis,
							fw_this_variant_this_seed_this_analysis,
							fw_this_variant_cohort_analysis,
							fw_variant_analysis]
						for compression in compression_fws:
							for data in data_fws:
								wf_links[data].append(compression)


				if BUILD_CAUSALITY_NETWORK:
					# BuildCausalityNetwork task
					fw_name = "BuildCausalityNetworkTask__Var_%d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)
					fw_this_variant_this_gen_this_sim_causality_network = Firework(
						BuildCausalityNetworkTask(
							input_results_directory = CELL_SIM_OUT_DIRECTORY,
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							output_dynamics_directory = CELL_SERIES_OUT_DIRECTORY,
							metadata = md_single,
							),
						name = fw_name,
						spec = {"_queueadapter": dict(analysis_q_cpus,
							job_name=fw_name), "_priority": 2}
						)

					wf_fws.append(fw_this_variant_this_gen_this_sim_causality_network)

					wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_this_gen_this_sim_causality_network)

					if COMPRESS_OUTPUT:
						# Don't compress any outputs or sim_data until
						# causality network scripts have finished running
						compression_fws = [
							fw_this_variant_sim_data_compression,
							fw_this_variant_this_gen_this_sim_compression
						]
						wf_links[fw_this_variant_this_gen_this_sim_causality_network].extend(compression_fws)

## Create workflow
log_info("Creating workflow.")

workflow = Workflow(wf_fws, links_dict = wf_links)

lpad.add_wf(workflow)
