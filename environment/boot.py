from __future__ import absolute_import, division, print_function

import errno
import os
import uuid
import argparse

import agent.event as event
from agent.outer import Outer
from agent.inner import Inner
from agent.shepherd import AgentShepherd
from agent.boot import EnvironmentControl

from environment.two_dim_lattice import EnvironmentSpatialLattice

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from models.ecoli.sim.simulation import EcoliSimulation

from wholecell.utils import constants
import wholecell.utils.filepath as fp
from models.ecoli.sim.variants import variant



class EnvironmentAgent(Outer):
	def build_state(self):
		lattice = {
			molecule: self.environment.lattice[index].tolist()
			for index, molecule in enumerate(self.environment.get_molecule_ids())}

		simulations = {
			agent_id: {
				'volume': simulation['state']['volume'],
				'location': self.environment.locations[agent_id][0:2].tolist(),
				'orientation': self.environment.locations[agent_id][2]}
			for agent_id, simulation in self.environment.simulations.iteritems()}

		return {
			'time': self.environment.time(),
			'lattice': lattice,
			'simulations': simulations}

	def update_state(self):
		self.send(self.kafka_config['environment_visualization'], self.build_state())

class BootEnvironmentSpatialLattice(object):
	def __init__(self, agent_id, agent_config):
		kafka_config = agent_config['kafka_config']
		raw_data = KnowledgeBaseEcoli()
		# create a dictionary with all saved environments
		self.environment_dict = {}
		for label in vars(raw_data.condition.environment):
			# initiate all molecules with 0 concentrations
			self.environment_dict[label] = {
				row["molecule id"]: 0 for row in raw_data.condition.environment_molecules
				}

			# get non-zero concentrations (assuming units.mmol / units.L)
			molecule_concentrations = getattr(raw_data.condition.environment, label)
			environment_non_zero_dict = {
				row["molecule id"]: row["concentration"].asNumber()
				for row in molecule_concentrations}

			# update environment_dict with non zero concentrations
			self.environment_dict[label].update(environment_non_zero_dict)

		# TODO (Eran) don't hardcode initial environment, get this from timeseries
		concentrations = self.environment_dict['minimal']

		self.environment = EnvironmentSpatialLattice(concentrations)
		self.outer = EnvironmentAgent(agent_id, kafka_config, self.environment)


class BootEcoli(object):
	'''
	This class initializes an EcoliSimulation, passes it to the `Inner` agent, and launches the simulation.
	The EcoliSimulation is initialized by passing it a pathname to sim_data, along with simulation parameters.
	'''
	def __init__(
			self, agent_id, agent_config,
			variant_type='wildtype', variant_index=0, seed=0):
		self.agent_id = agent_id
		kafka_config = agent_config['kafka_config']
		working_dir = agent_config['working_dir']
		outer_id = agent_config['outer_id']
		start_time = agent_config.get('start_time', 0)

		sim_path = fp.makedirs(working_dir, 'out', 'manual')
		sim_data_fit = os.path.join(sim_path, 'kb', 'simData_Most_Fit.cPickle')
		output_dir = os.path.join(sim_path, 'sim_' + self.agent_id, 'simOut')

		if not os.path.isfile(sim_data_fit):
			raise IOError(errno.ENOENT,
				'Missing "{}".  Run the Fitter?'.format(sim_data_fit))

		# Apply the variant to transform simData_Most_Fit.cPickle
		info, sim_data = variant.apply_variant(
			sim_data_file=sim_data_fit,
			variant_type=variant_type,
			variant_index=variant_index
			)

		options = {
			"simData":                sim_data,
			"outputDir":              output_dir,
			"initialTime":            start_time,
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
			"translationSupply":      True,
			}

		# Write a metadata file to aid analysis plots.
		# TODO(jerry): Skip it if another cell already wrote one?
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
			"translation_supply": options['translationSupply'],
			}
		metadata_dir = fp.makedirs(sim_path, 'metadata')
		metadata_path = os.path.join(metadata_dir, constants.JSON_METADATA_FILE)
		fp.write_json_file(metadata_path, metadata)

		self.simulation = EcoliSimulation(**options)
		self.inner = Inner(
			kafka_config,
			self.agent_id,
			outer_id,
			self.simulation)


class ShepherdControl(EnvironmentControl):

	"""
	Send messages to the other agents in the system to trigger execution and/or shutdown
	the Outer agent (which sends messages to shutdown all the associated Inner agents) or
	shutdown specific Inner agents directly (which then report back to the Outer agent and
	then terminate).
	"""

	def __init__(self, kafka_config):
		agent_id = 'shepherd_control'
		super(ShepherdControl, self).__init__(agent_id, kafka_config)

	# TODO (Ryan): set this up to send messages to a particular shepherd.
	def add_agent(self, agent_id, agent_type, agent_config):
		self.send(self.kafka_config['shepherd_control'], {
			'event': event.ADD_AGENT,
			'agent_id': agent_id,
			'agent_type': agent_type,
			'agent_config': agent_config})

	def remove_agent(self, prefix):
		""" Remove an agent given a prefix of its id """
		self.send(self.kafka_config['shepherd_control'], {
			'event': event.REMOVE_AGENT,
			'agent_prefix': prefix})

	def add_ecoli(self, lattice_id):
		self.add_agent(
			str(uuid.uuid1()),
			'ecoli',
			{'outer_id': lattice_id})

	def lattice_experiment(self, simulations):
		lattice_id = str(uuid.uuid1())
		self.add_agent(lattice_id, 'lattice', {})
		for index in range(simulations):
			self.add_ecoli(lattice_id)


def switch():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Run an agent for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=[
			'ecoli',
			'lattice',
			'shepherd',
			'experiment',
			'add',
			'remove',
			'trigger',
			'pause',
			'shutdown'],
		help='which command to boot')

	parser.add_argument(
		'--id',
		help='unique identifier for a new simulation agent')

	parser.add_argument(
		'--agent-id',
		help='unique identifier for addressing a particular agent')

	parser.add_argument(
		'--inner-id',
		help='unique identifier for inner agent')

	parser.add_argument(
		'--outer-id',
		help='unique identifier for outer agent')

	parser.add_argument(
		'-v', '--variant',
		default='wildtype',
		help='The variant type name. See models/ecoli/sim/variants/__init__.py'
			 ' for the choices')

	parser.add_argument(
		'-i', '--index',
		type=int,
		default=0,
		help='The variant index')

	parser.add_argument(
		'-s', '--seed',
		type=int,
		default=0,
		help='The simulation seed')

	parser.add_argument(
		'--kafka-host',
		default='127.0.0.1:9092',
		help='address for Kafka server')

	parser.add_argument(
		'--type',
		default='ecoli',
		help='type of agent to spawn in shepherd process')

	parser.add_argument(
		'--number',
		type=int,
		default=3,
		help='number of cell agents to spawn in lattice experiment')

	parser.add_argument(
		'--agent-receive',
		default='agent-receive',
		help='topic agents will receive control messages on')

	parser.add_argument(
		'--environment-control',
		default='environment-control',
		help='topic the environment will receive control messages on')

	parser.add_argument(
		'--simulation-receive',
		default='environment-broadcast',
		help='topic the simulations will receive messages on')

	parser.add_argument(
		'--simulation-send',
		default='environment-listen',
		help='topic the simulations will send messages on')

	parser.add_argument(
		'--environment-visualization',
		default='environment-state',
		help='topic the environment will send state information on')

	parser.add_argument(
		'--shepherd-control',
		default='shepherd-control',
		help='topic the shepherd will receive messages on')

	parser.add_argument(
		'--working-dir',
		default=os.getcwd(),
		help='the directory containing the sim path out/manual/'
	)

	args = parser.parse_args()
	kafka_config = {
		'host': args.kafka_host,
		'environment_control': args.environment_control,
		'simulation_receive': args.simulation_receive,
		'simulation_send': args.simulation_send,
		'environment_visualization': args.environment_visualization,
		'shepherd_control': args.shepherd_control,
		'subscribe_topics': []}

	if args.command == 'lattice':
		agent_id = args.id or 'lattice'
		BootEnvironmentSpatialLattice(agent_id, {'kafka_config': kafka_config})

	elif args.command == 'ecoli':
		if not args.id:
			raise ValueError('the "ecoli" command needs an --id argument')

		BootEcoli(
			args.id,
			{'kafka_config': kafka_config,
			 'working_dir': args.working_dir},
			variant_type=args.variant, variant_index=args.index, seed=args.seed)

	elif args.command == 'trigger':
		control = EnvironmentControl('environment_control', kafka_config)
		control.trigger_execution(args.id)
		control.shutdown()

	elif args.command == 'pause':
		control = EnvironmentControl('environment_control', kafka_config)
		control.pause_execution(args.id)
		control.shutdown()

	elif args.command == 'shutdown':
		control = EnvironmentControl('environment_control', kafka_config)
		control.shutdown_agent(args.id)
		control.shutdown()

	elif args.command == 'shepherd':
		initializers = {}

		def initialize_ecoli(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = kafka_config
			agent_config['working_dir'] = args.working_dir
			return BootEcoli(agent_id, agent_config)

		def initialize_lattice(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = kafka_config
			return BootEnvironmentSpatialLattice(agent_id, agent_config)

		initializers['lattice'] = initialize_lattice
		initializers['ecoli'] = initialize_ecoli

		shepherd = AgentShepherd('environment-shepherd', kafka_config, initializers)

	elif args.command == 'add':
		control = ShepherdControl(kafka_config)
		control.add_ecoli()
		control.shutdown()

	elif args.command == 'remove':
		control = ShepherdControl(kafka_config)
		control.remove_agent(args.id)
		control.shutdown()

	elif args.command == 'experiment':
		control = ShepherdControl(kafka_config)
		control.lattice_experiment(args.number)
		control.shutdown()

if __name__ == '__main__':
	switch()
