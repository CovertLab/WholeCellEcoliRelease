from __future__ import absolute_import, division, print_function

import os
import copy
import shutil

from lens.actor.inner import Inner
from lens.actor.boot import BootAgent
from lens.environment.make_media import Media
from lens.environment.boot import boot_lattice

from models.ecoli.sim.simulation import ecoli_simulation

from wholecell.utils import constants
import wholecell.utils.filepath as fp
from models.ecoli.sim.variants import apply_variant


def initialize_ecoli(config):
	'''
	Args:
		config (dict): options for initializing a simulation
	Returns:
		simulation (CellSimulation): The actual simulation which will perform the calculations.
	'''

	config = copy.deepcopy(config)

	outer_id = config['outer_id']
	sim_out_path = config['sim_out_path']
	variant_type = config['variant_type']
	variant_index = config['variant_index']

	sim_data_fit = os.path.join(sim_out_path, 'manual', 'kb', constants.SERIALIZED_SIM_DATA_FILENAME)
	variant_sim_data_directory = fp.makedirs(sim_out_path, 'agent', outer_id, 'kb')
	variant_sim_data_modified_file = os.path.join(variant_sim_data_directory, constants.SERIALIZED_SIM_DATA_MODIFIED)

	# copy sim_data into the experiment directory to support analysis
	# TODO (Eran) -- revisit this copy. Re-consider where to put the parca output.
	shutil.copy(sim_data_fit, variant_sim_data_modified_file)
	fp.verify_file_exists(sim_data_fit, 'Run runParca?')

	# Apply the variant to transform simData.cPickle
	info, sim_data_modified = apply_variant.apply_variant(
		sim_data_file=sim_data_fit,
		variant_type=variant_type,
		variant_index=variant_index)

	config['initialTime'] = config.pop('time') or 0
	config['simData'] = sim_data_modified
	return ecoli_simulation(**config)


def ecoli_boot_config(agent_config):
	'''
	Instantiates an initial or daughter EcoliSimulation, passes it to a new	`Inner` agent.
	Makes a simOut directory for the simulation in an embedded format:

		out/manual/experiment_id/cohort_id/generation_id/cell_id/simOut

	`agent_config` fields:
	    * generation (optional, the cell generation number)
	    * outer_id (id of outer environmental agent -- the experiment)
	    * working_dir (optional, wcEcoli path containing the sim path out/manual/)
	    * files (optional) list of data files:
			files[0] -- inherited_state_path to make a daughter cell
	    * start_time (optional)
	    * variant_type (optional)
	    * variant_index (optional)
	    * seed (optional)

	Returns:
		options (dict): simulation arguments for ecoli
	'''
	generation = agent_config.get('generation', 0)
	working_dir = agent_config.get('working_dir', os.getcwd())
	outer_id = agent_config.get('outer_id', 'lattice_000000')
	start_time = agent_config.get('start_time', 0)
	files = agent_config.get('files', [])
	inherited_state_path = files[0] if files else None
	variant_type = agent_config.get('variant_type', 'wildtype')
	variant_index = agent_config.get('variant_index', 0)
	seed = agent_config.get('seed', 0)
	volume = agent_config.get('volume', 1.0)
	index = agent_config.get('index', 0)

	# initialize state
	state = {
		'volume': volume,
		'environment_change': {}}
	agent_config['declare'] = state

	# TODO -- get actual cohort and cell ids
	# TODO -- change analysis scripts to allow the agent_id to be used here
	# TODO -- need to count number of initialized cells so that they won't over-write each other as 000000
	cohort_id = '%06d' % 0  # analysis scripts require starting with 0
	generation_id = 'generation_%06d' % generation
	cell_id = '%06d' % index # analysis scripts require starting with 0

	# make options for boot config
	sim_out_path = fp.makedirs(working_dir, 'out')
	output_dir = os.path.join(sim_out_path, 'agent', outer_id, cohort_id, generation_id, cell_id, 'simOut')
	metadata_dir = fp.makedirs(sim_out_path, 'agent', 'metadata')
	metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)

	options = {
		"sim_out_path":           sim_out_path,
		"variant_type":           variant_type,
		"variant_index":          variant_index,
		"outer_id":               outer_id,
		"lengthSec":              3*60*60,
		"outputDir":              output_dir,
		"initialTime":            start_time,
		"inheritedStatePath":     inherited_state_path,
		"logToDisk":              True,
		"overwriteExistingFiles": True,
		"seed":                   seed,
		"timeStepSafetyFraction": 1.3,
		"maxTimeStep":            0.9,
		"updateTimeStepFreq":     5,
		"logToShell":             True,
		"logToDiskEvery":         1,
		"massDistribution":       True,
		"growthRateNoise":        False,
		"dPeriodDivision":        False,
		"translationSupply":      True}

	# Write a metadata file to aid analysis plots.
	metadata = {
		"git_hash":           fp.run_cmdline("git rev-parse HEAD"),
		"git_branch":         fp.run_cmdline("git symbolic-ref --short HEAD"),
		"description":        "an Ecoli Cell Agent",
		"time":               fp.timestamp(),
		"total_gens":         0,  # not known in advance for multi-scale sims
		"analysis_type":      None,
		"variant":            variant_type,
		"mass_distribution":  options['massDistribution'],
		"growth_rate_noise":  options['growthRateNoise'],
		"d_period_division":  options['dPeriodDivision'],
		"translation_supply": options['translationSupply']}
	fp.write_json_file(metadata_path, metadata)

	return options

def boot_ecoli(agent_id, agent_type, agent_config):
	'''
	Pass configuration and initialization function to a new `Inner` agent, which will launch
	the simulation.
	'''

	if 'outer_id' not in agent_config:
		raise ValueError("--outer-id required")

	agent_config['boot_config'] = ecoli_boot_config(agent_config)

	inner = Inner(
		agent_id,
		agent_type,
		agent_config,
		initialize_ecoli)

	return inner



class BootEcoli(BootAgent):
	def __init__(self):
		super(BootEcoli, self).__init__()
		self.agent_types = {
			'lattice': boot_lattice,
			'ecoli': boot_ecoli}

def run():
	boot = BootEcoli()
	boot.execute()

if __name__ == '__main__':
	run()
