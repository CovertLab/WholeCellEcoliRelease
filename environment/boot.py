from __future__ import absolute_import, division, print_function

import os
import uuid
import cPickle
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

from wholecell.fireworks.firetasks import VariantSimDataTask



class EnvironmentAgent(Outer):
	def build_state(self):
		lattice = {
			molecule: self.environment.lattice[index].tolist()
			for index, molecule in enumerate(self.environment.get_molecule_ids())}

		simulations = {
			agent_id: {
				'volume': state['volume'],
				'location': self.environment.locations[agent_id][0:2].tolist(),
				'orientation': self.environment.locations[agent_id][2]}
			for agent_id, state in self.environment.simulations.iteritems()}

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
	def __init__(self, agent_id, agent_config):
		self.agent_id = agent_id
		kafka_config = agent_config['kafka_config']
		working_dir = agent_config['working_dir']

		sim_data_fit = os.path.join(os.getcwd(),'out','manual','kb','simData_Most_Fit.cPickle')
		sim_data_variant = os.path.join(os.getcwd(), 'out', 'manual', 'kb', 'simData_Modified.cPickle')
		variant_metadata = os.path.join(os.getcwd(), 'out', 'manual', 'metadata')
		output_dir = os.path.join(os.getcwd(), 'out', 'manual', 'sim_' + self.agent_id,'simOut')

		# copy the file simData_Most_Fit.cPickle to simData_Modified.cPickle
		task = VariantSimDataTask(
			variant_function='wildtype',
			variant_index=0,
			input_sim_data=sim_data_fit,
			output_sim_data=sim_data_variant,
			variant_metadata_directory=variant_metadata,
		)
		task.run_task({})

		with open(sim_data_variant, "rb") as input_sim_data:
			sim_data = cPickle.load(input_sim_data)

		options = {}
		options["simData"] = sim_data
		options["outputDir"] = output_dir
		options["logToDisk"] = True
		options["overwriteExistingFiles"] = True

		options["seed"] = 0
		options["lengthSec"] = 10800
		options["timeStepSafetyFraction"] = 1.3
		options["maxTimeStep"] = 0.9
		options["updateTimeStepFreq"] = 5
		options["logToShell"] = True
		options["logToDiskEvery"] = 1
		options["massDistribution"] = True
		options["growthRateNoise"] = False
		options["dPeriodDivision"] = False
		options["translationSupply"] = True

		self.simulation = EcoliSimulation(**options)
		self.inner = Inner(
			kafka_config,
			self.agent_id,
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

	def add_agent(self, agent_id, agent_type, agent_config):
		self.send(self.kafka_config['shepherd_control'], {
			'event': event.ADD_AGENT,
			'agent_id': agent_id,
			'agent_type': agent_type,
			'agent_config': agent_config})

	def add_ecoli(self):
		self.add_agent(str(uuid.uuid1()), 'ecoli', {})

	def remove_agent(self, prefix):
		""" Remove an agent given a prefix of its id """
		self.send(self.kafka_config['shepherd_control'], {
			'event': event.REMOVE_AGENT,
			'agent_prefix': prefix})

	def lattice_experiment(self, simulations):
		self.add_agent('lattice', 'lattice', {})
		for index in range(simulations):
			self.add_ecoli()


def switch():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	# TODO (eran) share argparse code with agent/boot.py
	# One way to do that is via a base class similar to ScriptBase.py or shared subroutines.
	# another way is argparse.ArgumentParser(parents=[parent_parser]): https://docs.python.org/2/library/argparse.html?highlight=argparse#parents

	parser = argparse.ArgumentParser(
		description='Run the agents for the environmental context simulation')

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
			'shutdown'],
		help='which command to boot')

	parser.add_argument(
		'--id',
		help='unique identifier for simulation agent')

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
		help='number of agents to spawn in lattice experiment')

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
		help='the directory containing the project files'
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
		BootEnvironmentSpatialLattice('lattice', {'kafka_config': kafka_config})

	elif args.command == 'ecoli':
		if not args.id:
			raise ValueError('--id must be supplied for ecoli command')

		BootEcoli(args.id, {
			'kafka_config': kafka_config,
			'working_dir': args.working_dir})

	elif args.command == 'trigger':
		control = EnvironmentControl('environment_control', kafka_config)
		control.trigger_execution()
		control.shutdown()

	elif args.command == 'shutdown':
		control = EnvironmentControl('environment_control', kafka_config)

		if not args.id:
			control.shutdown_environment()
		else:
			control.shutdown_simulation(args.id)
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

		shepherd = AgentShepherd('shepherd', kafka_config, initializers)

	elif args.command == 'add':
		control = ShepherdControl(kafka_config)
		control.add_ecoli()

	elif args.command == 'remove':
		control = ShepherdControl(kafka_config)
		control.remove_agent(args.id)

	elif args.command == 'experiment':
		control = ShepherdControl(kafka_config)
		control.lattice_experiment(args.number)

if __name__ == '__main__':
	switch()
