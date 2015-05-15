#!/usr/bin/env python

from fireworks import Firework, LaunchPad, Workflow, ScriptTask
from wholecell.fireworks.firetasks import InitKbTask
from wholecell.fireworks.firetasks import SymlinkTask
from wholecell.fireworks.firetasks import FitKbTask
from wholecell.fireworks.firetasks import VariantKbTask
from wholecell.fireworks.firetasks import SimulationTask
from wholecell.fireworks.firetasks import SimulationDaughterTask
from wholecell.fireworks.firetasks import AnalysisSingleTask
from wholecell.fireworks.firetasks import AnalysisMultiGenTask
from wholecell.sim.simulation import DEFAULT_SIMULATION_KWARGS
from jinja2 import Template

import wholecell.utils.constants
import yaml
import os
import datetime
import subprocess
import collections

def run_cmd(cmd):
	environ = {
	"PATH": os.environ["PATH"],
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

### Set path variables

dirname = os.path.dirname
WC_ECOLI_DIRECTORY = dirname(dirname(os.path.abspath(__file__)))
OUT_DIRECTORY = os.path.join(WC_ECOLI_DIRECTORY, "out")

now = datetime.datetime.now()
SUBMISSION_TIME = "%04d%02d%02d.%02d%02d%02d.%06d" % (
	now.year, now.month, now.day,
	now.hour, now.minute, now.second,
	now.microsecond)
KB_DIRECTORY = os.path.join(OUT_DIRECTORY, SUBMISSION_TIME, "kb")
METADATA_DIRECTORY = os.path.join(OUT_DIRECTORY, SUBMISSION_TIME, "metadata")

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
N_INIT_SIMS = int(os.environ.get("N_INIT_SIMS", "1"))
N_GENS = int(os.environ.get("N_GENS", "1"))
SINGLE_DAUGHTERS = bool(int(os.environ.get("SINGLE_DAUGHTERS", "0")))

### Create directories

if not os.path.exists(OUT_DIRECTORY):
	os.makedirs(OUT_DIRECTORY)

if not os.path.exists(KB_DIRECTORY):
	os.makedirs(KB_DIRECTORY)

if not os.path.exists(METADATA_DIRECTORY):
	os.makedirs(METADATA_DIRECTORY)

for i in VARIANTS_TO_RUN:
	VARIANT_DIRECTORY = os.path.join(OUT_DIRECTORY, SUBMISSION_TIME, VARIANT + "_%06d" % i)
	VARIANT_KB_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")

	if not os.path.exists(VARIANT_DIRECTORY):
		os.makedirs(VARIANT_DIRECTORY)

	if not os.path.exists(VARIANT_KB_DIRECTORY):
		os.makedirs(VARIANT_KB_DIRECTORY)

	if not os.path.exists(VARIANT_METADATA_DIRECTORY):
		os.makedirs(VARIANT_METADATA_DIRECTORY)

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

write_file(
	os.path.join(METADATA_DIRECTORY, "git_hash"),
	run_cmd(["git", "rev-parse", "HEAD"])
	)

write_file(
	os.path.join(METADATA_DIRECTORY, "git_branch"),
	run_cmd(["git", "symbolic-ref", "--short", "HEAD"])
	)

write_file(
	os.path.join(METADATA_DIRECTORY, "git_diff"),
	run_cmd(["git", "diff"])
	)

write_file(
	os.path.join(METADATA_DIRECTORY, "description"),
	os.environ.get("DESC", "")
	)


#### Create workflow

# Create launchpad
lpad = LaunchPad(**yaml.load(open("my_launchpad.yaml")))

# Store list of FireWorks
wf_fws = []

# Store links defining parent/child relationships of FireWorks
wf_links = collections.defaultdict(list)


### Initialize KB

filename_kb_fit_0 = (
			wholecell.utils.constants.SERIALIZED_KB_PREFIX +
			"_Fit_0" +
			wholecell.utils.constants.SERIALIZED_KB_SUFFIX
			)

fw_name = "InitKbTask"
fw_initKb = Firework(
	InitKbTask(
		output = os.path.join(KB_DIRECTORY, filename_kb_fit_0)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_initKb)

# Unfit KB compression

fw_name = "ScriptTask_compression_fit_0_KB"
fw_kb_fit_0_compression = Firework(
	ScriptTask(
		script = "bzip2 " + os.path.join(KB_DIRECTORY, filename_kb_fit_0)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_kb_fit_0_compression)

## Create symlink to unfit KB

fw_name = "SymlinkTask_KB_Unfit"
fw_symlink_unfit = Firework(
	SymlinkTask(
		to = filename_kb_fit_0,
		link = os.path.join(KB_DIRECTORY, wholecell.utils.constants.SERIALIZED_KB_UNFIT_FILENAME),
		overwrite_if_exists = True
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_symlink_unfit)

wf_links[fw_initKb].append(fw_symlink_unfit)


### Fit (Level 1)

filename_kb_fit_1 = (
			wholecell.utils.constants.SERIALIZED_KB_PREFIX +
			"_Fit_1" +
			wholecell.utils.constants.SERIALIZED_KB_SUFFIX
			)

fw_name = "FitKbTask_Level_1"
fw_fit_level_1 = Firework(
	FitKbTask(
		fit_level = 1,
		input_kb = os.path.join(KB_DIRECTORY, filename_kb_fit_0),
		output_kb = os.path.join(KB_DIRECTORY, filename_kb_fit_1),
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_fit_level_1)
wf_links[fw_symlink_unfit].append(fw_fit_level_1)

# Fit Level 1 KB compression

fw_name = "ScriptTask_compression_fit_1_KB"
fw_kb_fit_1_compression = Firework(
	ScriptTask(
		script = "bzip2 " + os.path.join(KB_DIRECTORY, filename_kb_fit_1)
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_kb_fit_1_compression)

## Create symlink to most fit KB
# (when more fitting stages are implemented, move this down)

fw_name = "SymlinkTask_KB_Most_Fit"
fw_symlink_most_fit = Firework(
	SymlinkTask(
		to = filename_kb_fit_1,
		link = os.path.join(KB_DIRECTORY, wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME),
		overwrite_if_exists = True
		),
	name = fw_name,
	spec = {"_queueadapter": {"job_name": fw_name}}
	)

wf_fws.append(fw_symlink_most_fit)

wf_links[fw_fit_level_1].append(fw_symlink_most_fit)


### Create variants and simulations
for i in VARIANTS_TO_RUN:
	VARIANT_DIRECTORY = os.path.join(OUT_DIRECTORY, SUBMISSION_TIME, VARIANT + "_%06d" % i)
	VARIANT_KB_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "kb")
	VARIANT_METADATA_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "metadata")

	# Variant KB creation task
	fw_name = "VariantKbTask"
	fw_this_variant_kb = Firework(
		VariantKbTask(
			variant_function = VARIANT,
			variant_index = i,
			input_kb = os.path.join(KB_DIRECTORY, wholecell.utils.constants.SERIALIZED_KB_MOST_FIT_FILENAME),
			output_kb = os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle"),
			variant_metadata_directory = VARIANT_METADATA_DIRECTORY,
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}}
		)

	wf_fws.append(fw_this_variant_kb)

	wf_links[fw_symlink_most_fit].append(fw_this_variant_kb)

	# Variant KB compression
	fw_name = "ScriptTask_compression_variant_KB"
	fw_this_variant_kb_compression = Firework(
		ScriptTask(
			script = "bzip2 " + os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle")
			),
		name = fw_name,
		spec = {"_queueadapter": {"job_name": fw_name}}
		)

	wf_fws.append(fw_this_variant_kb_compression)

	for j in xrange(N_INIT_SIMS):
		SEED_DIRECTORY = os.path.join(VARIANT_DIRECTORY, "%06d" % j)
		SEED_PLOT_DIRECTORY = os.path.join(SEED_DIRECTORY, "plotOut")

		fw_name = "AnalysisMultiGenTask__Seed_%06d" % (j)
		fw_this_variant_this_seed_this_analysis = Firework(
			AnalysisMultiGenTask(
				input_seed_directory = SEED_DIRECTORY,
				input_kb = os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle"),
				output_plots_directory = SEED_PLOT_DIRECTORY,
				),
			name = fw_name,
			spec = {"_queueadapter": {"job_name": fw_name}}
			)
		wf_fws.append(fw_this_variant_this_seed_this_analysis)

		wf_links[fw_this_variant_this_seed_this_analysis].append(fw_this_variant_kb_compression)
		wf_links[fw_this_variant_this_seed_this_analysis].append(fw_kb_fit_0_compression) # Maybe not necessary
		wf_links[fw_this_variant_this_seed_this_analysis].append(fw_kb_fit_1_compression) # Maybe not necessary

		sims_this_seed = collections.defaultdict(list)

		for k in xrange(N_GENS):
			GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % k)

			for l in (xrange(2**k) if not SINGLE_DAUGHTERS else [0]):
				CELL_DIRECTORY = os.path.join(GEN_DIRECTORY, "%06d" % l)
				CELL_SIM_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "simOut")
				CELL_PLOT_OUT_DIRECTORY = os.path.join(CELL_DIRECTORY, "plotOut")

				# TODO: Add conditional logic here for mother vs daughter cells
				# Simulation task
				fw_name = "SimulationTask__Gen_%d__Cell_%d" % (k, l)

				if k == 0:
					fw_this_variant_this_gen_this_sim = Firework(
						SimulationTask(
							input_kb = os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle"),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							seed = j,
							length_sec = WC_LENGTHSEC,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name}}
						)
				elif k > 0:
					PARENT_GEN_DIRECTORY = os.path.join(SEED_DIRECTORY, "generation_%06d" % (k - 1))
					PARENT_CELL_DIRECTORY = os.path.join(PARENT_GEN_DIRECTORY, "%06d" % (l // 2))
					PARENT_CELL_SIM_OUT_DIRECTORY = os.path.join(PARENT_CELL_DIRECTORY, "simOut")
					DAUGHTER_STATE_DIRECTORY = os.path.join(PARENT_CELL_SIM_OUT_DIRECTORY, "Daughter%d" % (l % 2 + 1))

					fw_this_variant_this_gen_this_sim = Firework(
						SimulationDaughterTask(
							input_kb = os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle"),
							output_directory = CELL_SIM_OUT_DIRECTORY,
							inherited_state_path = DAUGHTER_STATE_DIRECTORY,
							seed = (j + 1) * ((2**k - 1) + l),
							length_sec = WC_LENGTHSEC,
							),
						name = fw_name,
						spec = {"_queueadapter": {"job_name": fw_name}}
						)

				wf_fws.append(fw_this_variant_this_gen_this_sim)
				wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_this_seed_this_analysis)

				sims_this_seed[k].append(fw_this_variant_this_gen_this_sim)

				if k == 0:
					wf_links[fw_this_variant_kb].append(fw_this_variant_this_gen_this_sim)

				elif k > 0:
					fw_parent_sim = sims_this_seed[k - 1][l // 2]
					wf_links[fw_parent_sim].append(fw_this_variant_this_gen_this_sim)

				# AnalysisSingle task
				fw_name = "AnalysisSingleTask__Gen_%d__Cell_%d" % (k, l)
				fw_this_variant_this_gen_this_sim_analysis = Firework(
					AnalysisSingleTask(
						input_results_directory = CELL_SIM_OUT_DIRECTORY,
						input_kb = os.path.join(VARIANT_KB_DIRECTORY, "KnowledgeBase_Modified.cPickle"),
						output_plots_directory = CELL_PLOT_OUT_DIRECTORY,
						),
					name = fw_name,
					spec = {"_queueadapter": {"job_name": fw_name}}
					)

				wf_fws.append(fw_this_variant_this_gen_this_sim_analysis)

				wf_links[fw_this_variant_this_gen_this_sim].append(fw_this_variant_this_gen_this_sim_analysis)

				# NOTE: Change the following line if you are debugging and want to skip analysis
				# For example, make it:
				# last_fw_that_needs_kbs = fw_this_variant_this_gen_this_sim
				last_fw_that_needs_kbs = fw_this_variant_this_gen_this_sim_analysis

				wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_this_variant_kb_compression)
				wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_kb_fit_0_compression) # Maybe not necessary
				wf_links[fw_this_variant_this_gen_this_sim_analysis].append(fw_kb_fit_1_compression) # Maybe not necessary


### Create workflow

workflow = Workflow(wf_fws, links_dict = wf_links)

lpad.add_wf(workflow)
