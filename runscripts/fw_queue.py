#!/usr/bin/env python

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

import wholecell.utils.constants
import yaml
import os
import datetime
import subprocess
import collections
import cPickle

def run_cmd(cmd):
	environ = {
		"PATH": os.environ["PATH"],
		"LD_LIBRARY_PATH": os.environ["LD_LIBRARY_PATH"],
		"LANG": "C",
		"LC_ALL": "C",
		}
	out = subprocess.Popen(cmd, stdout = subprocess.PIPE, env=environ).communicate()[0]
	return out

def write_file(filename, content):
	h = open(filename, "w")
	h.write(content)
	h.close()

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

### Set other environment variables

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

if not RUN_AGGREGATE_ANALYSIS:
	COMPRESS_OUTPUT = False

### Set path variables

dirname = os.path.dirname
WC_ECOLI_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
OUT_DIRECTORY = os.path.join(WC_ECOLI_DIRECTORY, "out")
CACHED_SIM_DATA_DIRECTORY = os.path.join(WC_ECOLI_DIRECTORY, "cached")

now = datetime.datetime.now()
SUBMISSION_TIME = "%04d%02d%02d.%02d%02d%02d.%06d" % (
	now.year, now.month, now.day,
	now.hour, now.minute, now.second,
	now.microsecond)
INDIV_OUT_DIRECTORY = os.path.join(OUT_DIRECTORY, SUBMISSION_TIME + "__" + SIM_DESCRIPTION)
KB_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, "kb")
METADATA_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, "metadata")

### Create directories

if not os.path.exists(OUT_DIRECTORY):
	os.makedirs(OUT_DIRECTORY)

if not os.path.exists(INDIV_OUT_DIRECTORY):
	os.makedirs(INDIV_OUT_DIRECTORY)

if not os.path.exists(KB_DIRECTORY):
	os.makedirs(KB_DIRECTORY)

if not os.path.exists(METADATA_DIRECTORY):
	os.makedirs(METADATA_DIRECTORY)

if VERBOSE_QUEUE:
	print "Building filestructure."

for i in VARIANTS_TO_RUN:
	VARIANT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")
	VARIANT_COHORT_PLOT_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "plotOut")

	if not os.path.exists(VARIANT_DIRECTORY):
		os.makedirs(VARIANT_DIRECTORY)

	if not os.path.exists(VARIANT_SIM_DATA_DIRECTORY):
		os.makedirs(VARIANT_SIM_DATA_DIRECTORY)

	if not os.path.exists(VARIANT_METADATA_DIRECTORY):
		os.makedirs(VARIANT_METADATA_DIRECTORY)

	if not os.path.exists(VARIANT_COHORT_PLOT_DIRECTORY):
		os.makedirs(VARIANT_COHORT_PLOT_DIRECTORY)

	for j in xrange(N_INIT_SIMS):
		SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, "plotOut")

		if not os.path.exists(SEED_DIRECTORY):
			os.makedirs(SEED_DIRECTORY)

		if not os.path.exists(SEED_PLOT_DIRECTORY):
			os.makedirs(SEED_PLOT_DIRECTORY)

		for k in xrange(N_GENS):
			GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)

			if not os.path.exists(GEN_DIRECTORY):
				os.makedirs(GEN_DIRECTORY)

			for l in (xrange(2**k) if not SINGLE_DAUGHTERS else [0]):
				CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "plotOut")

				if not os.path.exists(CELL_DIRECTORY):
					os.makedirs(CELL_DIRECTORY)

				if not os.path.exists(CELL_SIM_OUT_DIRECTORY):
					os.makedirs(CELL_SIM_OUT_DIRECTORY)

				if not os.path.exists(CELL_PLOT_OUT_DIRECTORY):
					os.makedirs(CELL_PLOT_OUT_DIRECTORY)

### Write metadata
metadata = {
	"git_hash": run_cmd(["git", "rev-parse", "HEAD"]),
	"git_branch": run_cmd(["git", "symbolic-ref", "--short", "HEAD"]),
	"git_diff": run_cmd(["git", "diff"]),
	"description": os.environ.get("DESC", ""),
	"time": SUBMISSION_TIME,
	"total_gens": str(N_GENS),
	"analysis_type": None,
	"variant": VARIANT,
	"mass_distribution" : MASS_DISTRIBUTION,
	"growth_rate_noise" : GROWTH_RATE_NOISE,
	"d_period_division" : D_PERIOD_DIVISION,
	"translation_supply" : TRANSLATION_SUPPLY,
	}

for key, value in metadata.iteritems():
	if type(value) != str:
		continue
	write_file(os.path.join(METADATA_DIRECTORY, key), value)

h = open(os.path.join(METADATA_DIRECTORY, "metadata.cPickle"), "w")
cPickle.dump(metadata, h, cPickle.HIGHEST_PROTOCOL)
h.close()

#### Create workflow

# Create launchpad
lpad = LaunchPad(**yaml.load(open(LAUNCHPAD_FILE)))

# Store list of FireWorks
wf_fws = []

# Store links defining parent/child relationships of FireWorks
wf_links = collections.defaultdict(list)


### Initialize KB

filename_raw_data = wholecell.utils.constants.SERIALIZED_RAW_DATA

fw_name = "InitRawData"

if VERBOSE_QUEUE:
	print "Queuing {}".format(fw_name)

fw_init_raw_data = Firework(
	InitRawDataTask(
		output = os.path.join(KB_DIRECTORY, filename_raw_data)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_init_raw_data)

### Fit (Level 1)

filename_sim_data_fit_1 = (
			wholecell.utils.constants.SERIALIZED_SIM_DATA_PREFIX +
			"_Fit_1" +
			wholecell.utils.constants.SERIALIZED_SIM_DATA_SUFFIX
			)

fw_name = "FitSimDataTask_Level_1"

if VERBOSE_QUEUE:
	print "Queuing {}".format(fw_name)

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
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name, "cpus_per_task": cpusForFitter}, "_priority":1}
	)

wf_fws.append(fw_fit_level_1)
wf_links[fw_init_raw_data].append(fw_fit_level_1)

# Unfit KB compression
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_raw_data"

	if VERBOSE_QUEUE:
		print "Queuing {}".format(fw_name)

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

if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_sim_data_1"

	if VERBOSE_QUEUE:
		print "Queuing {}".format(fw_name)

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
	print "Queuing {}".format(fw_name)

fw_symlink_most_fit = Firework(
	SymlinkTask(
		to = filename_sim_data_fit_1,
		link = os.path.join(KB_DIRECTORY, wholecell.utils.constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME),
		overwrite_if_exists = True
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
	)

wf_fws.append(fw_symlink_most_fit)

wf_links[fw_fit_level_1].append(fw_symlink_most_fit)


### Initialize validation data

# Initiate raw validation data
filename_raw_validation_data = wholecell.utils.constants.SERIALIZED_RAW_VALIDATION_DATA

fw_name = "InitValidationDataRaw"

if VERBOSE_QUEUE:
	print "Queuing {}".format(fw_name)

fw_raw_validation_data = Firework(
	InitRawValidationDataTask(
		output = os.path.join(KB_DIRECTORY, filename_raw_validation_data)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}, "_priority":1}
	)

wf_fws.append(fw_raw_validation_data)

# Raw validation data compression
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
filename_validation_data = (
			wholecell.utils.constants.SERIALIZED_VALIDATION_DATA
			)

fw_name = "InitValidationData"

if VERBOSE_QUEUE:
	print "Queuing {}".format(fw_name)

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

# Full validation data compression
if COMPRESS_OUTPUT:
	fw_name = "ScriptTask_compression_validation_data"

	if VERBOSE_QUEUE:
		print "Queuing {}".format(fw_name)

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
for i in VARIANTS_TO_RUN:
	if VERBOSE_QUEUE:
		print "Queuing Variant {}".format(i)
	VARIANT_DIRECTORY = os.path.join(INDIV_OUT_DIRECTORY, VARIANT + "_%06d" % i)
	VARIANT_SIM_DATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")
	metadata["variant_function"] = VARIANT
	metadata["variant_index"] = i

	# Variant simData creation task
	fw_name = "VariantSimDataTask_%06d" % (i)
	fw_this_variant_sim_data = Firework(
		VariantSimDataTask(
			variant_function = VARIANT,
			variant_index = i,
			input_sim_data = os.path.join(KB_DIRECTORY, wholecell.utils.constants.SERIALIZED_SIM_DATA_MOST_FIT_FILENAME),
			output_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
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
				script = "bzip2 -v " + os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle")
				),
			name = fw_name,
			spec = {"_queueadapter": {"job_name": fw_name}, "_priority":0}
			)

		wf_fws.append(fw_this_variant_sim_data_compression)

	# Cohort analysis
	COHORT_PLOT_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "plotOut")

	if RUN_AGGREGATE_ANALYSIS:
		metadata["analysis_type"] = "cohort"
		fw_name = "AnalysisCohortTask__Var_%02d" % (i)
		fw_this_variant_cohort_analysis = Firework(
			AnalysisCohortTask(
				input_variant_directory = VARIANT_DIRECTORY,
				input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
				input_validation_data = os.path.join(KB_DIRECTORY, filename_validation_data),
				output_plots_directory = COHORT_PLOT_DIRECTORY,
				metadata = metadata,
				),
			name = fw_name,
			spec = {"_queueadapter": {"job_name": fw_name}, "_priority":4}
			)
		wf_fws.append(fw_this_variant_cohort_analysis)

	for j in xrange(N_INIT_SIMS):
		if VERBOSE_QUEUE:
			print "\tQueuing Seed {}".format(j)
		SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, "plotOut")
		metadata["seed"] = j

		if RUN_AGGREGATE_ANALYSIS:
			metadata["analysis_type"] = 'multigen'

			fw_name = "AnalysisMultiGenTask__Var_%02d__Seed_%06d" % (i, j)
			fw_this_variant_this_seed_this_analysis = Firework(
				AnalysisMultiGenTask(
					input_seed_directory = SEED_DIRECTORY,
					input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
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
				print "\t\tQueuing Gen %02d." % (k)
			GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)
			metadata["gen"] = k

			for l in (xrange(2**k) if not SINGLE_DAUGHTERS else [0]):

				if VERBOSE_QUEUE:
					print "\t\t\tQueuing Cell {}".format(l)
				CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "plotOut")

				# TODO: Add conditional logic here for mother vs daughter cells
				# Simulation task
				fw_name = "SimulationTask__Var_%02d__Seed_%d__Gen_%d__Cell_%d" % (i, j, k, l)

				if k == 0:
					fw_this_variant_this_gen_this_sim = Firework(
						SimulationTask(
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
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
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
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

				wf_fws.append(fw_this_variant_this_gen_this_sim)
				if RUN_AGGREGATE_ANALYSIS:
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
							input_sim_data = os.path.join(VARIANT_SIM_DATA_DIRECTORY, "simData_Modified.cPickle"),
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
						wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_this_variant_sim_data_compression)
						wf_links[fw_this_variant_this_seed_this_analysis].append(fw_this_variant_sim_data_compression)
						wf_links[fw_this_variant_cohort_analysis].append(fw_this_variant_sim_data_compression)
						wf_links[fw_variant_analysis].append(fw_this_variant_sim_data_compression)
						wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_validation_data_compression)
						wf_links[fw_this_variant_this_seed_this_analysis].append(fw_validation_data_compression)
						wf_links[fw_this_variant_cohort_analysis].append(fw_validation_data_compression)
						wf_links[fw_variant_analysis].append(fw_validation_data_compression)
						wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_this_variant_this_gen_this_sim_compression)
						wf_links[fw_this_variant_this_seed_this_analysis].append(fw_this_variant_this_gen_this_sim_compression)
						wf_links[fw_this_variant_cohort_analysis].append(fw_this_variant_this_gen_this_sim_compression)
						wf_links[fw_variant_analysis].append(fw_this_variant_this_gen_this_sim_compression)

## Create workflow


if VERBOSE_QUEUE:
	print "Creating workflow."

workflow = Workflow(wf_fws, links_dict = wf_links)

lpad.add_wf(workflow)
