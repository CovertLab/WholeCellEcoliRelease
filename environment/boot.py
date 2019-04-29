from __future__ import absolute_import, division, print_function

import os
import errno

from agent.outer import Outer
from agent.inner import Inner
from agent.boot import BootAgent

from environment.lattice import EnvironmentSpatialLattice
from environment.surrogates.chemotaxis import Chemotaxis
from environment.surrogates.endocrine import Endocrine
from environment.surrogates.transport_lookup_minimal import TransportMinimal
from models.ecoli.sim.simulation import ecoli_simulation
from environment.condition.make_media import Media

from wholecell.utils import constants
import wholecell.utils.filepath as fp
from models.ecoli.sim.variants import apply_variant

DEFAULT_COLOR = [0.6, 0.4, 0.3]

class EnvironmentAgent(Outer):
	def build_state(self):
		lattice = {
			molecule: self.environment.lattice[index].tolist()
			for index, molecule in enumerate(self.environment.get_molecule_ids())}

		simulations = {
			agent_id: {
				'volume': simulation['state']['volume'],
				'color': simulation['state'].get('color', DEFAULT_COLOR),
				'location': self.environment.locations[agent_id][0:2].tolist(),
				'orientation': self.environment.locations[agent_id][2],
				'parent_id': simulation.get('parent_id', '')}
			for agent_id, simulation in self.environment.simulations.iteritems()}

		return {
			'outer_id': self.agent_id,
			'agent_type': 'ecoli',
			'running': not self.paused,
			'time': self.environment.time(),
			'edge_length': self.environment.edge_length,
			'cell_radius': self.environment.cell_radius,
			'lattice': lattice,
			'simulations': simulations}

	def update_state(self):
		self.send(
			self.topics['visualization_receive'],
			self.build_state(),
			print_send=False)

def boot_lattice(agent_id, agent_type, agent_config):
	media_id = agent_config.get('media_id', 'minimal')
	media = agent_config.get('media', {})
	print("Media condition: {}".format(media_id))
	if not media:
		make_media = Media()
		media = make_media.make_recipe(media_id)

	agent_config['concentrations'] = media
	environment = EnvironmentSpatialLattice(agent_config)

	return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

# wcEcoli initialize and boot
def initialize_ecoli(boot_config, synchronize_config):
	'''
	Args:
		boot_config (dict): options for initializing a simulation
		synchronize_config (dict): additional options that can be passed in for initialization
	Returns:
		simulation (CellSimulation): The actual simulation which will perform the calculations.
	'''
	synchronize_config['initialTime'] = synchronize_config.pop('time')
	boot_config.update(synchronize_config)
	return ecoli_simulation(**boot_config)

def get_ecoli_boot_config(agent_id, agent_config):
	'''
	Sets up an ecoli simulation. Makes the metadata files. Returns a dictionary of options
	Args:
		agent_id (str): the id of this agent
		agent_config (dict) with fields
			* kafka_config
			* outer_id (id of environmental context agent)
			* working_dir (optional, wcEcoli path containing the sim path out/manual/)
			* files (optional) list of data files:
				files[0] -- inherited_state_path to make a daughter cell
			* start_time (optional)
			* variant_type (optional)
			* variant_index (optional)
			* seed (optional)
			* volume (optional)
	Returns:
		options (dict): simulation arguments for ecoli
	'''

	working_dir = agent_config.get('working_dir', os.getcwd())
	outer_id = agent_config['outer_id']
	start_time = agent_config.get('start_time', 0)
	files = agent_config.get('files', [])
	inherited_state_path = files[0] if files else None
	variant_type = agent_config.get('variant_type', 'wildtype')
	variant_index = agent_config.get('variant_index', 0)
	seed = agent_config.get('seed', 0)

	# initialize state
	state = {
		'volume': 1.0,
		'environment_change': {}}
	agent_config['state'] = state

	# make options for boot config
	sim_path = fp.makedirs(working_dir, 'out', 'manual')
	sim_data_fit = os.path.join(sim_path, 'kb', 'simData_Most_Fit.cPickle')
	output_dir = os.path.join(sim_path, 'sim_' + agent_id, 'simOut')

	if not os.path.isfile(sim_data_fit):
		raise IOError(
			errno.ENOENT,
			'Missing "{}".  Run the Parca?'.format(sim_data_fit))

	# Apply the variant to transform simData_Most_Fit.cPickle
	info, sim_data = apply_variant.apply_variant(
		sim_data_file=sim_data_fit,
		variant_type=variant_type,
		variant_index=variant_index)

	options = {
		"simData":                sim_data,
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
		# "total_gens":       1,  # not known in advance for multi-scale sims
		"analysis_type":      None,
		"variant":            variant_type,
		"mass_distribution":  options['massDistribution'],
		"growth_rate_noise":  options['growthRateNoise'],
		"d_period_division":  options['dPeriodDivision'],
		"translation_supply": options['translationSupply']}
	metadata_dir = fp.makedirs(sim_path, 'metadata')
	metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
	fp.write_json_file(metadata_path, metadata)

	return options


def boot_ecoli(agent_id, agent_type, agent_config):
	'''
	passes configuration and initialization function to a new `Inner` agent, which will launch the simulation.
	'''
	if 'outer_id' not in agent_config:
		raise ValueError("--outer-id required")
	outer_id = agent_config['outer_id']
	options = get_ecoli_boot_config(agent_id, agent_config)

	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		options,
		initialize_ecoli)

	return inner

# Chemotaxis surrogate initialize and boot
def initialize_chemotaxis(boot_config, synchronize_config):
	'''
	Args:
		boot_config (dict): options for initializing a simulation
		synchronize_config (dict): additional options that can be passed in for initialization
	Returns:
		simulation (CellSimulation): The actual simulation which will perform the calculations.
	'''
	boot_config.update(synchronize_config)
	return Chemotaxis(boot_config)

def boot_chemotaxis(agent_id, agent_type, agent_config):
	agent_id = agent_id
	outer_id = agent_config['outer_id']

	# initialize state and options
	state = {
		'volume': 1.0,
		'environment_change': {}}
	agent_config['state'] = state
	options = {}

	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		options,
		initialize_chemotaxis)

	return inner

# Endocrine surrogate initialize and boot
def initialize_endocrine(boot_config, synchronize_config):
	'''
	Args:
		boot_config (dict): options for initializing a simulation
		synchronize_config (dict): additional options that can be passed in for initialization
	Returns:
		simulation (CellSimulation): The actual simulation which will perform the calculations.
	'''
	boot_config.update(synchronize_config)
	return Endocrine(boot_config)

def boot_endocrine(agent_id, agent_type, agent_config):
	agent_id = agent_id
	outer_id = agent_config['outer_id']

	# initialize state and options
	state = {
		'volume': 1.0,
		'environment_change': {}}
	agent_config['state'] = state
	options = {}

	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		options,
		initialize_endocrine)

	return inner

# Transport lookup minimal surrogate initialize and boot
def initialize_transport_minimal(boot_config, synchronize_config):
	boot_config.update(synchronize_config)
	return TransportMinimal(boot_config)

def boot_transport_minimal(agent_id, agent_type, agent_config):
	agent_id = agent_id
	outer_id = agent_config['outer_id']

	# initialize state and options
	state = {
		'volume': 1.0,
		'environment_change': {}}
	agent_config['state'] = state
	options = {}

	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		options,
		initialize_transport_minimal)

	return inner


class BootEnvironment(BootAgent):
	def __init__(self):
		super(BootEnvironment, self).__init__()
		self.agent_types = {
			'lattice': boot_lattice,
			'ecoli': boot_ecoli,
			'chemotaxis': boot_chemotaxis,
			'endocrine': boot_endocrine,
			'transport_minimal': boot_transport_minimal,
			}

if __name__ == '__main__':
	boot = BootEnvironment()
	boot.execute()
