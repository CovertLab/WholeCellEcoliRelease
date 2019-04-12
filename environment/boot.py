from __future__ import absolute_import, division, print_function

import os
import time
import errno

import agent.event as event
from agent.outer import Outer
from agent.inner import Inner
from agent.boot import BootAgent

from environment.lattice import EnvironmentSpatialLattice
from environment.surrogates.chemotaxis import Chemotaxis

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

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
	media = agent_config.get('media', 'minimal')
	print("Media condition: {}".format(media))
	raw_data = KnowledgeBaseEcoli()

	# make media object
	make_media = Media()
	new_media = make_media.make_recipe(media)

	agent_config['concentrations'] = new_media
	environment = EnvironmentSpatialLattice(agent_config)

	return EnvironmentAgent(agent_id, agent_type, agent_config, environment)

def boot_ecoli(agent_id, agent_type, agent_config):
	'''
	Instantiates an initial or daughter EcoliSimulation, passes it to a new
	`Inner` agent, and launches the simulation. `agent_config` fields:
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
	'''
	if 'outer_id' not in agent_config:
		raise ValueError("--outer-id required")

	kafka_config = agent_config['kafka_config']
	working_dir = agent_config.get('working_dir', os.getcwd())
	outer_id = agent_config['outer_id']
	start_time = agent_config.get('start_time', 0)
	files = agent_config.get('files', [])
	inherited_state_path = files[0] if files else None
	variant_type = agent_config.get('variant_type', 'wildtype')
	variant_index = agent_config.get('variant_index', 0)
	seed = agent_config.get('seed', 0)
	volume = agent_config.get('volume', 1.2)

	# Create the inner agent before instantiating the simulation so it can send
	# a message to the lattice without waiting for the simulation to initialize.
	# TODO(jerry): Is the agent OK receiving messages w/o a sim? Add a
	#    setter for its simulation property, make it queue messages until it
	#    gets one, change the constructor type signature to allow None?
	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		None)

	inner.send(kafka_config['topics']['environment_receive'], {
		'event': event.CELL_DECLARE,
		'agent_id': outer_id,
		'inner_id': agent_id,
		'agent_config': agent_config,
		'state': {
			'volume': volume,
			'environment_change': {}}})

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

	inner.simulation = ecoli_simulation(**options)

	return inner

def boot_chemotaxis(agent_id, agent_type, agent_config):
	agent_id = agent_id
	outer_id = agent_config['outer_id']
	volume = 1.0
	kafka_config = agent_config['kafka_config']

	inner = Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		None)

	inner.send(kafka_config['topics']['environment_receive'], {
		'event': event.CELL_DECLARE,
		'agent_id': outer_id,
		'inner_id': agent_id,
		'agent_config': agent_config,
		'state': {
			'volume': volume,
			'environment_change': {}}})

	simulation = Chemotaxis()
	inner.simulation = simulation

	time.sleep(5)  # TODO(jerry): Wait for the Chemotaxis to boot

	return inner


class BootEnvironment(BootAgent):
	def __init__(self):
		super(BootEnvironment, self).__init__()
		self.agent_types = {
			'lattice': boot_lattice,
			'ecoli': boot_ecoli,
			'chemotaxis': boot_chemotaxis,
			}

if __name__ == '__main__':
	boot = BootEnvironment()
	boot.execute()
