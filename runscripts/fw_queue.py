#!/usr/bin/env python

'''
Creates an array of firetasks (wf_fws) and specifies their links (wf_links) in
a workflow (wf) for Fireworks, and submits them to the queue.

Several environmental variables can be specified, shown below with their (type, default value).
These are set as follows, and otherwise revert to their default value:

	VARIABLE=VALUE python runscripts/fw_queue.py

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
		to run the simulation instead of running the fitter; useful for repeated
		simulations where raw_data/fitter files are not changed
	PARALLEL_FITTER (int, "0"): if nonzero, some fitter operations will run in
		parallel instead of serially
	DEBUG_FITTER (int, "0"): if nonzero, this reduces the number of TFs and
		conditions; allows for faster debugging of fitter
	COMPRESS_OUTPUT (int, "0"): if nonzero, outputs will be compressed (.bz2)
	RUN_AGGREGATE_ANALYSIS (int, "1"): if nonzero, all analyses are run on
		simulation output
	DISABLE_RIBOSOME_CAPACITY_FITTING (int, "0"): if nonzero, ribosome
		expression is not fit to protein synthesis demands
	DISABLE_RNAPOLY_CAPACITY_FITTING (int, "0"): if nonzero, RNA polymerase
		expression is not fit to RNA synthesis demands
	ADJUST_RNA_AND_PROTEIN_PARAMETERS (int, "1"): if nonzero, some RNA and
		protein expression parameters will be adjusted to get expression
	ADJUST_RNASE_EXPRESSION (int, "0"): if nonzero, adjusts the expression of all
		RNase mRNA lower
	DISABLE_MEASURED_PROTEIN_DEG (int, "0"): if nonzero, does not use any measured
		protein degradation rates and defaults to the N-end rule

Simulation parameters:
	N_GENS (int, "1"): the number of generations to be simulated
	N_INIT_SIMS (int, "1"): the number of initial simulations
	SINGLE_DAUGHTERS (int, "0"): if nonzero, the simulation will generate only
		one daughter cell for each new generation rather than two, thus avoiding
		an exponential increase in the number of simulations
	WC_LENGTHSEC (int, "10800"): sets the maximum simulation time in seconds, useful
		for short simulations (default is 3 hr)
	TIMESTEP_MAX (float, "0.9"): sets the maximum time step
	TIMESTEP_SAFETY_FRAC (float, "1.3"): increases the time step by this factor
		if conditions are favorable; up the the limit of the max time step
	TIMESTEP_UPDATE_FREQ (int, "5"): frequency at which the time step is updated

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

Additional variables:
	LAUNCHPAD_FILE (str, "my_launchpad.yaml"): set launchpad config file location
	VERBOSE_QUEUE (int, "1"): if nonzero, gives more detailed messages during
		fireworks set up

Environment variables that matter when running the workflow:
	DEBUG_GC (int, "0"): if nonzero, enable leak detection in the analysis plots
	WC_ANALYZE_FAST (anything, --): if set, run each analysis plot in a separate
		process
'''

from fireworks import Firework, LaunchPad, Workflow, ScriptTask
from wholecell.fireworks.firetasks import InitRawDataTask
from wholecell.fireworks.firetasks import InitRawValidationDataTask
from wholecell.fireworks.firetasks import InitValidationDataTask
from wholecell.fireworks.firetasks import SymlinkTask
from wholecell.fireworks.firetasks import FitSimDataTask
from wholecell.fireworks.firetasks import VariantSimDataTask
from wholecell.fireworks.firetasks import SimulationTask
from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.fireworks.firetasks import AnalysisVariantTask
from wholecell.fireworks.firetasks import AnalysisCohortTask
from wholecell.fireworks.firetasks import AnalysisSingleTask
from wholecell.fireworks.firetasks import AnalysisMultiGenTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS

from wholecell.utils import constants
from wholecell.utils import filepath
import yaml
import os
import collections
import cPickle


#### Initial setup ###


### Set variant variables

VARIANT = os.environ.get("VARIANT", "wildtype")
FIRST_VARIANT_INDEX = int(os.environ.get("FIRST_VARIANT_INDEX", "0"))
LAST_VARIANT_INDEX = int(os.environ.get("LAST_VARIANT_INDEX", "0"))

if LAST_VARIANT_INDEX == -1:
	from models.ecoli.sim.variants import nameToNumIndicesMapping
	LAST_VARIANT_INDEX = nameToNumIndicesMapping[VARIANT]

# This variable gets iterated over in multiple places
# So be careful if you change it to xrange
VARIANTS_TO_RUN = range(FIRST_VARIANT_INDEX, LAST_VARIANT_INDEX + 1)

### Set other simulation parameters

WC_LENGTHSEC = int(os.environ.get("WC_LENGTHSEC", DEFAULT_SIMULATION_KWARGS["lengthSec"]))
TIMESTEP_SAFETY_FRAC = float(os.environ.get("TIMESTEP_SAFETY_FRAC", DEFAULT_SIMULATION_KWARGS["timeStepSafetyFraction"]))
TIMESTEP_MAX = float(os.environ.get("TIMESTEP_MAX", DEFAULT_SIMULATION_KWARGS["maxTimeStep"]))
TIMESTEP_UPDATE_FREQ = int(os.environ.get("TIMESTEP_UPDATE_FREQ", DEFAULT_SIMULATION_KWARGS["updateTimeStepFreq"]))
MASS_DISTRIBUTION = bool(int(os.environ.get("MASS_DISTRIBUTION", DEFAULT_SIMULATION_KWARGS["massDistribution"])))
GROWTH_RATE_NOISE = bool(int(os.environ.get("GROWTH_RATE_NOISE", DEFAULT_SIMULATION_KWARGS["growthRateNoise"])))
D_PERIOD_DIVISION = bool(int(os.environ.get("D_PERIOD_DIVISION", DEFAULT_SIMULATION_KWARGS["dPeriodDivision"])))
TRANSLATION_SUPPLY = bool(int(os.environ.get("TRANSLATION_SUPPLY", DEFAULT_SIMULATION_KWARGS["translationSupply"])))
N_INIT_SIMS = int(os.environ.get("N_INIT_SIMS", "1"))
N_GENS = int(os.environ.get("N_GENS", "1"))
SINGLE_DAUGHTERS = bool(int(os.environ.get("SINGLE_DAUGHTERS", "0")))
LAUNCHPAD_FILE = str(os.environ.get("LAUNCHPAD_FILE", "my_launchpad.yaml"))
COMPRESS_OUTPUT = bool(int(os.environ.get("COMPRESS_OUTPUT", "0")))
SIM_DESCRIPTION = os.environ.get("DESC", "").replace(" ", "_")
VERBOSE_QUEUE = bool(int(os.environ.get("VERBOSE_QUEUE", "1")))
RUN_AGGREGATE_ANALYSIS = bool(int(os.environ.get("RUN_AGGREGATE_ANALYSIS", "1")))
CACHED_SIM_DATA = bool(int(os.environ.get("CACHED_SIM_DATA", "0")))
PARALLEL_FITTER = bool(int(os.environ.get("PARALLEL_FITTER", "0")))
DEBUG_FITTER = bool(int(os.environ.get("DEBUG_FITTER", "0")))
DISABLE_RIBOSOME_CAPACITY_FITTING = bool(int(os.environ.get("DISABLE_RIBOSOME_CAPACITY_FITTING", "0")))
DISABLE_RNAPOLY_CAPACITY_FITTING = bool(int(os.environ.get("DISABLE_RNAPOLY_CAPACITY_FITTING", "0")))
RNAPOLY_ACTIVITY_FITTING = bool(int(os.environ.get("RNAPOLY_ACTIVITY_FITTING", "0")))
MRNA_HALF_LIFE_FITTING = bool(int(os.environ.get("MRNA_HALF_LIFE_FITTING", "0")))
MAX_RNAP_ACTIVITY = bool(int(os.environ.get("MAX_RNAP_ACTIVITY", "0")))
ADJUST_RNA_AND_PROTEIN_PARAMETERS = bool(int(os.environ.get("ADJUST_RNA_AND_PROTEIN_PARAMETERS", "1")))
ADJUST_RNASE_EXPRESSION = bool(int(os.environ.get("ADJUST_RNASE_EXPRESSION", "0")))
DISABLE_MEASURED_PROTEIN_DEG = bool(int(os.environ.get("DISABLE_MEASURED_PROTEIN_DEG", "0")))
VARIABLE_ELONGATION_TRANSCRIPTION = bool(int(os.environ.get("VARIABLE_ELONGATION_TRANSCRIPTION", "0")))
VARIABLE_ELONGATION_TRANSLATION = bool(int(os.environ.get("VARIABLE_ELONGATION_TRANSLATION", "0")))
ALTERNATE_MASS_FRACTION_PROTEIN = bool(int(os.environ.get("ALTERNATE_MASS_FRACTION_PROTEIN", "0")))
ALTERNATE_MASS_FRACTION_RNA = bool(int(os.environ.get("ALTERNATE_MASS_FRACTION_RNA", "0")))
ALTERNATE_MASS_FRACTION_MRNA = bool(int(os.environ.get("ALTERNATE_MASS_FRACTION_MRNA", "0")))
ALTERNATE_R_PROTEIN_DEGRADATION = bool(int(os.environ.get("ALTERNATE_R_PROTEIN_DEGRADATION", "0")))
ALTERNATE_RNA_SEQ = bool(int(os.environ.get("ALTERNATE_RNA_SEQ", "0")))
ALTERNATE_RNA_HALF_LIFE = bool(int(os.environ.get("ALTERNATE_RNA_HALF_LIFE", "0")))
ALTERNATE_TRANSLATION_EFFICIENCY = bool(int(os.environ.get("ALTERNATE_TRANSLATION_EFFICIENCY", "0")))
ALTERNATE_RIBOSOME_ACTIVITY = bool(int(os.environ.get("ALTERNATE_RIBOSOME_ACTIVITY", "0")))
DISABLE_RNAP_FRACTION_INCREASE = bool(int(os.environ.get("DISABLE_RNAP_FRACTION_INCREASE", "0")))
DISABLE_RIBOSOME_ACTIVITY_FIX = bool(int(os.environ.get("DISABLE_RIBOSOME_ACTIVITY_FIX", "0")))
SAVE_CELL_SPECS = bool(int(os.environ.get("SAVE_CELL_SPECS", "0")))
CELL_SPECS_FILE = bool(int(os.environ.get("CELL_SPECS_FILE", "0")))
WRITE_TRANSLATION_EFFICIENCIES = bool(int(os.environ.get("WRITE_TRANSLATION_EFFICIENCIES", "0")))

if not RUN_AGGREGATE_ANALYSIS:
	COMPRESS_OUTPUT = False

### Set path variables and create directories

WC_ECOLI_DIRECTORY = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT_DIRECTORY = filepath.makedirs(WC_ECOLI_DIRECTORY, "out")
CACHED_SIM_DATA_DIRECTORY = os.path.join(WC_ECOLI_DIRECTORY, "cached")

SUBMISSION_TIME = filepath.timestamp()
INDIV_OUT_DIRECTORY = filepath.makedirs(OUT_DIRECTORY, SUBMISSION_TIME + "__" + SIM_DESCRIPTION)
KB_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, "kb")
METADATA_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, "metadata")


if VERBOSE_QUEUE:
	print "Building filestructure."

for i in VARIANTS_TO_RUN:
	VARIANT_DIRECTORY = filepath.makedirs(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "metadata")
	VARIANT_COHORT_PLOT_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "plotOut")

	for j in xrange(N_INIT_SIMS):
		SEED_DIRECTORY = filepath.makedirs(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, "plotOut")

		for k in xrange(N_GENS):
			GEN_DIRECTORY = filepath.makedirs(SEED_DIRECTORY, "generation_%06d" % k)

			for l in (xrange(2**k) if not SINGLE_DAUGHTERS else [0]):
				CELL_DIRECTORY = filepath.makedirs(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = filepath.makedirs(CELL_DIRECTORY, "plotOut")


### Write metadata
metadata = {
	"git_hash": filepath.run_cmdline("git rev-parse HEAD"),
	"git_branch": filepath.run_cmdline("git symbolic-ref --short HEAD"),
	"git_diff": filepath.run_cmdline("git diff", trim=False),
	"description": os.environ.get("DESC", ""),
	"time": SUBMISSION_TIME,
	"total_gens": str(N_GENS),
	"analysis_type": None,
	"variant": VARIANT,
	"mass_distribution": MASS_DISTRIBUTION,
	"growth_rate_noise": GROWTH_RATE_NOISE,
	"d_period_division": D_PERIOD_DIVISION,
	"translation_supply": TRANSLATION_SUPPLY,
	}

for key, value in metadata.iteritems():
	if not isinstance(value, basestring):
		continue
	filepath.write_file(os.path.join(METADATA_DIRECTORY, key), value)

with open(os.path.join(METADATA_DIRECTORY, constants.SERIALIZED_METADATA_FILE), "wb") as f:
	cPickle.dump(metadata, f, cPickle.HIGHEST_PROTOCOL)

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
filename_sim_data_modified = constants.SERIALIZED_SIM_DATA_MODIFIED

fw_name = "InitRawData"

if VERBOSE_QUEUE:
	print "Queueing {}".format(fw_name)

fw_init_raw_data = Firework(
	InitRawDataTask(
		output = os.path.join(KB_DIRECTORY, filename_raw_data)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_init_raw_data)

### Fit (Level 1)

filename_sim_data_fit_1 = constants.SERIALIZED_FIT1_FILENAME

fw_name = "FitSimDataTask_Level_1"

if VERBOSE_QUEUE:
	print "Queueing {}".format(fw_name)

if PARALLEL_FITTER:
	cpusForFitter = 8
else:
	cpusForFitter = 1
fw_fit_level_1 = Firework(
	FitSimDataTask(
		fit_level = 1,
		input_data = os.path.join(KB_DIRECTORY, filename_raw_data),
		output_data = os.path.join(KB_DIRECTORY, filename_sim_data_fit_1),
		cached = CACHED_SIM_DATA,
		cached_data = os.path.join(CACHED_SIM_DATA_DIRECTORY, filename_sim_data_fit_1),
		cpus = cpusForFitter,
		debug = DEBUG_FITTER,
		disable_ribosome_capacity_fitting = DISABLE_RIBOSOME_CAPACITY_FITTING,
		disable_rnapoly_capacity_fitting = DISABLE_RNAPOLY_CAPACITY_FITTING,
		rnapoly_activity_fitting = RNAPOLY_ACTIVITY_FITTING,
		mrna_half_life_fitting = MRNA_HALF_LIFE_FITTING,
		max_rnap_activity = MAX_RNAP_ACTIVITY,
		variable_elongation_transcription = VARIABLE_ELONGATION_TRANSCRIPTION,
		variable_elongation_translation = VARIABLE_ELONGATION_TRANSLATION,
		adjust_rna_and_protein_parameters = ADJUST_RNA_AND_PROTEIN_PARAMETERS,
		adjust_rnase_expression = ADJUST_RNASE_EXPRESSION,
		disable_measured_protein_deg = DISABLE_MEASURED_PROTEIN_DEG,
		alternate_mass_fraction_protein = ALTERNATE_MASS_FRACTION_PROTEIN,
		alternate_mass_fraction_rna = ALTERNATE_MASS_FRACTION_RNA,
		alternate_mass_fraction_mrna = ALTERNATE_MASS_FRACTION_MRNA,
		alternate_r_protein_degradation = ALTERNATE_R_PROTEIN_DEGRADATION,
		alternate_rna_seq = ALTERNATE_RNA_SEQ,
		alternate_rna_half_life = ALTERNATE_RNA_HALF_LIFE,
		alternate_translation_efficiency = ALTERNATE_TRANSLATION_EFFICIENCY,
		alternate_ribosome_activity = ALTERNATE_RIBOSOME_ACTIVITY,
		disable_rnap_fraction_increase = DISABLE_RNAP_FRACTION_INCREASE,
		disable_ribosome_activity_fix = DISABLE_RIBOSOME_ACTIVITY_FIX,
		save_cell_specs = SAVE_CELL_SPECS,
		cell_specs_file = CELL_SPECS_FILE,
		write_translation_efficiencies = WRITE_TRANSLATION_EFFICIENCIES
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": cpusForFitter}, "_priority":1}
	)

wf_fws.append(fw_fit_level_1)
wf_links[fw_init_raw_data].append(fw_fit_level_1)

# Unfit KB compression
fw_raw_data_compression = None
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_raw_data"

	if VERBOSE_QUEUE:
		print "Queueing {}".format(fw_name)

	fw_raw_data_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_raw_data)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_raw_data_compression)
	wf_links[fw_fit_level_1].append(fw_raw_data_compression)

# Fit Level 1 KB compression

fw_sim_data_1_compression = None
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_sim_data_1"

	if VERBOSE_QUEUE:
		print "Queueing {}".format(fw_name)

	fw_sim_data_1_compression = Firework(
		ScriptTask(
			script = "bzip2 -v " + os.path.join(KB_DIRECTORY, filename_sim_data_fit_1)
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
		)

	wf_fws.append(fw_sim_data_1_compression)

## Create symlink to most fit KB
# (when more fitting stages are implemented, move this down)

fw_name = "SymlinkTask_KB_Most_Fit"

if VERBOSE_QUEUE:
	print "Queueing {}".format(fw_name)

fw_symlink_most_fit = Firework(
	SymlinkTask(
		to = filename_sim_data_fit_1,
		link = os.path.join(KB_DIRECTORY, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME),
		overwrite_if_exists = True
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
	)

wf_fws.append(fw_symlink_most_fit)

wf_links[fw_fit_level_1].append(fw_symlink_most_fit)


### Initialize validation data

# Initiate raw validation data
filename_raw_validation_data = constants.SERIALIZED_RAW_VALIDATION_DATA

fw_name = "InitValidationDataRaw"

if VERBOSE_QUEUE:
	print "Queueing {}".format(fw_name)

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

if VERBOSE_QUEUE:
	print "Queueing {}".format(fw_name)

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

	if VERBOSE_QUEUE:
		print "Queueing {}".format(fw_name)

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

if RUN_AGGREGATE_ANALYSIS:

	metadata["analysis_type"] = "variant"
	metadata["total_variants"] = str(len(VARIANTS_TO_RUN))

	fw_name = "AnalysisVariantTask"
	fw_variant_analysis = Firework(
		AnalysisVariantTask(
			input_directory = os.path.join(INDIV_OUT_DIRECTORY),
			input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
			output_plots_directory = VARIANT_PLOT_DIRECTORY,
			metadata = metadata,
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":5}
		)
	wf_fws.append(fw_variant_analysis)

### Create variants and simulations
fw_this_variant_sim_data_compression = None
fw_this_variant_this_gen_this_sim_compression = None

for i in VARIANTS_TO_RUN:
	if VERBOSE_QUEUE:
		print "Queueing Variant {} {}".format(VARIANT, i)
	VARIANT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")
	metadata["variant_function"] = VARIANT
	metadata["variant_index"] = i

	# Variant simData creation task
	fw_name = "VariantSimDataTask_%06d" % (i,)
	fw_this_variant_sim_data = Firework(
		VariantSimDataTask(
			variant_function = VARIANT,
			variant_index = i,
			input_sim_data = os.path.join(KB_DIRECTORY, constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME),
			output_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
			variant_metadata_directory = VARIANT_METADATA_DIRECTORY,
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
		)

	wf_fws.append(fw_this_variant_sim_data)

	wf_links[fw_symlink_most_fit].append(fw_this_variant_sim_data)

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
		metadata["analysis_type"] = "cohort"
		fw_name = "AnalysisCohortTask__Var_%02d" % (i,)
		fw_this_variant_cohort_analysis = Firework(
			AnalysisCohortTask(
				input_variant_directory = VARIANT_DIRECTORY,
				input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
				input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
				output_plots_directory = COHORT_PLOT_DIRECTORY,
				metadata = metadata,
				),
			name = fw_name,
			spec = {"_queueadapter": {"job_name": fw_name}, "_priority":4}
			)
		wf_fws.append(fw_this_variant_cohort_analysis)

	fw_this_variant_this_seed_this_analysis = None

	for j in xrange(N_INIT_SIMS):
		if VERBOSE_QUEUE:
			print "\tQueueing Seed {}".format(j)
		SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, "plotOut")
		metadata["seed"] = j

		if RUN_AGGREGATE_ANALYSIS:
			metadata["analysis_type"] = 'multigen'

			fw_name = "AnalysisMultiGenTask__Var_%02d__Seed_%06d" % (i, j)
			fw_this_variant_this_seed_this_analysis = Firework(
				AnalysisMultiGenTask(
					input_seed_directory = SEED_DIRECTORY,
					input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
					input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
					output_plots_directory = SEED_PLOT_DIRECTORY,
					metadata = metadata,
					),
				name = fw_name,
				spec = {"_queueadapter": {"job_name": fw_name}, "_priority":3}
				)
			wf_fws.append(fw_this_variant_this_seed_this_analysis)

			if COMPRESS_OUTPUT:
				wf_links[fw_this_variant_this_seed_this_analysis].append(fw_this_variant_sim_data_compression)

		sims_this_seed = collections.defaultdict(list)

		for k in xrange(N_GENS):
			if VERBOSE_QUEUE:
				print "\t\tQueueing Gen %02d." % (k,)
			GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)
			metadata["gen"] = k

			for l in (xrange(2**k) if not SINGLE_DAUGHTERS else [0]):

				if VERBOSE_QUEUE:
					print "\t\t\tQueueing Cell {}".format(l)
				CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "plotOut")

				# TODO: Add conditional logic here for mother vs daughter cells
				# Simulation task
				fw_name = "SimulationTask__Var_%02d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)

				if k == 0:
					fw_this_variant_this_gen_this_sim = Firework(
						SimulationTask(
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							seed = j,
							length_sec = WC_LENGTHSEC,
							timestep_safety_frac = TIMESTEP_SAFETY_FRAC,
							timestep_max = TIMESTEP_MAX,
							timestep_update_freq = TIMESTEP_UPDATE_FREQ,
							mass_distribution = MASS_DISTRIBUTION,
							growth_rate_noise = GROWTH_RATE_NOISE,
							d_period_division = D_PERIOD_DIVISION,
							translation_supply = TRANSLATION_SUPPLY,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": 1}, "_priority":10}
						)
				elif k > 0:
					PARENT_GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % (k - 1))
					PARENT_CELL_DIRECTORY = os.path.join(PARENT_GEN_DIRECTORY, "%06d" % (l // 2))
					PARENT_CELL_SIM_OUT_DIRECTORY = os.path.join(PARENT_CELL_DIRECTORY, "simOut")
					DAUGHTER_STATE_DIRECTORY = os.path.join(PARENT_CELL_SIM_OUT_DIRECTORY, "Daughter%d" % (l % 2 + 1))

					fw_this_variant_this_gen_this_sim = Firework(
						SimulationDaughterTask(
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							inherited_state_path = DAUGHTER_STATE_DIRECTORY,
							seed = (j + 1) * ((2**k - 1) + l),
							length_sec = WC_LENGTHSEC,
							timestep_safety_frac = TIMESTEP_SAFETY_FRAC,
							timestep_max = TIMESTEP_MAX,
							timestep_update_freq = TIMESTEP_UPDATE_FREQ,
							mass_distribution = MASS_DISTRIBUTION,
							growth_rate_noise = GROWTH_RATE_NOISE,
							d_period_division = D_PERIOD_DIVISION,
							translation_supply = TRANSLATION_SUPPLY,
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

				if RUN_AGGREGATE_ANALYSIS:
					if COMPRESS_OUTPUT:
						# Output compression job
						fw_name = "ScriptTask_compression_simulation__Seed_%d__Gen_%d__Cell_%d" % (j, k, l)
						fw_this_variant_this_gen_this_sim_compression = Firework(
							ScriptTask(
								script = 'for dir in %s; do echo "Compressing $dir"; find "$dir" -type f | xargs bzip2; done' % os.path.join(CELL_SIM_OUT_DIRECTORY, "*")
								),
							name = fw_name,
							spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
							)

						wf_fws.append(fw_this_variant_this_gen_this_sim_compression)

					metadata["analysis_type"] = "single"

					# AnalysisSingle task
					fw_name = "AnalysisSingleTask__Var_%d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)
					fw_this_variant_this_gen_this_sim_analysis = Firework(
						AnalysisSingleTask(
							input_results_directory = CELL_SIM_OUT_DIRECTORY,
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, filename_sim_data_modified),
							input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
							output_plots_directory = CELL_PLOT_OUT_DIRECTORY,
							metadata = metadata,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name}, "_priority":2}
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

## Create workflow


if VERBOSE_QUEUE:
	print "Creating workflow."

workflow = Workflow(wf_fws, links_dict = wf_links)

lpad.add_wf(workflow)
