from __future__ import absolute_import, division, print_function

import errno
import os
import uuid
import argparse

import agent.event as event
from agent.outer import Outer
from agent.inner import Inner
from agent.shepherd import AgentShepherd
from agent.boot import EnvironmentControl, AgentCommand

from environment.two_dim_lattice import EnvironmentSpatialLattice

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from models.ecoli.sim.simulation import ecoli_simulation

from wholecell.utils import constants
import wholecell.utils.filepath as fp
from models.ecoli.sim.variants import apply_variant



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
	Instantiates an initial or daughter EcoliSimulation, passes it to the
	`Inner` agent, and launches the simulation. `agent_config` fields:
		* kafka_config
		* working_dir (wcEcoli path containing the sim path out/manual/)
		* inherited_state_path (optional, to make a daughter cell)
		* variant_type (optional)
		* variant_index (optional)
		* seed (optional)
	'''
	def __init__(
			self, agent_id, agent_config,
			variant_type='wildtype', variant_index=0, seed=0):
		self.agent_id = agent_id

		kafka_config = agent_config['kafka_config']
		working_dir = agent_config['working_dir']
		outer_id = agent_config['outer_id']
		start_time = agent_config.get('start_time', 0)
		inherited_state_path = agent_config.get('inherited_state_path', None)
		variant_type = agent_config.get('variant_type', 'wildtype')
		variant_index = agent_config.get('variant_index', 0)
		seed = agent_config.get('seed', 0)

		sim_path = fp.makedirs(working_dir, 'out', 'manual')
		sim_data_fit = os.path.join(sim_path, 'kb', 'simData_Most_Fit.cPickle')
		output_dir = os.path.join(sim_path, 'sim_' + self.agent_id, 'simOut')

		if not os.path.isfile(sim_data_fit):
			raise IOError(errno.ENOENT,
				'Missing "{}".  Run the Fitter?'.format(sim_data_fit))

		# Apply the variant to transform simData_Most_Fit.cPickle
		info, sim_data = apply_variant.apply_variant(
			sim_data_file=sim_data_fit,
			variant_type=variant_type,
			variant_index=variant_index
			)

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

		self.simulation = ecoli_simulation(**options)
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

	def add_ecoli(self, agent_config):
		self.add_agent(
			str(uuid.uuid1()),
			'ecoli',
			agent_config)

	def lattice_experiment(self, simulations):
		lattice_id = str(uuid.uuid1())
		self.add_agent(lattice_id, 'lattice', {})
		for index in range(simulations):
			self.add_ecoli({'outer_id': lattice_id})


class EnvironmentCommand(AgentCommand):
	"""
	Extend `AgentCommand` with new commands related to the lattice and ecoli experiments
	"""

	def __init__(self):
		choices = ['ecoli', 'lattice']
		description = '''
		Run an agent for the environmental context simulation.
		The commands are:
		`add [--variant V] [--index I] [--seed S]` ask the Shepherd to add an Ecoli agent,
		`ecoli --id ID [--working-dir D] [--variant V] [--index I] [--seed S]` run an Ecoli agent in this process,
		`experiment [--number N]` ask the Shepherd to run a Lattice agent and N Ecoli agents,
		`lattice` run a Lattice environment agent in this process,
		`pause` pause the simulation,
		`remove --id ID` ask all Shepherds to remove agents "ID*",
		`shepherd [--working-dir D]` run a Shepherd agent in this process,
		`shutdown` shut down the environment agent and all its connected agents,
		`trigger` start running the simulation'''

		super(EnvironmentCommand, self).__init__(choices, description)

	def add_arguments(self, parser):
		parser = super(EnvironmentCommand, self).add_arguments(parser)

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

		return parser

	def shepherd_initializers(self, args):
		initializers = super(EnvironmentCommand, self).shepherd_initializers(args)

		def initialize_ecoli(agent_id, agent_config):
			agent_config = dict(agent_config,
				kafka_config=self.kafka_config,
				working_dir=args.working_dir)
			return BootEcoli(agent_id, agent_config)

		def initialize_lattice(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = self.kafka_config
			return BootEnvironmentSpatialLattice(agent_id, agent_config)

		initializers['lattice'] = initialize_lattice
		initializers['ecoli'] = initialize_ecoli

		return initializers

	def lattice(self, args):
		agent_id = args.id or 'lattice'
		BootEnvironmentSpatialLattice(agent_id, {'kafka_config': self.kafka_config})

	def ecoli(self, args):
		if not args.id:
			raise ValueError('the "ecoli" command needs an --id argument')
		if not args.outer_id:
			raise ValueError('the "ecoli" command needs an --outer-id argument')

		agent_config = dict(
			kafka_config=self.kafka_config,
			working_dir=args.working_dir,
			variant_type=args.variant,
			variant_index=args.index,
			seed=args.seed,
			outer_id=args.outer_id
			)
		BootEcoli(args.id, agent_config)

	def add(self, args):
		agent_config = dict(
			variant_type=args.variant,
			variant_index=args.index,
			seed=args.seed,
			outer_id=args.id,
			)
		control = ShepherdControl(self.kafka_config)
		control.add_ecoli(agent_config)
		control.shutdown()

	def remove(self, args):
		if not args.id:
			raise ValueError('the "remove" command needs an --id argument')
		control = ShepherdControl(self.kafka_config)
		control.remove_agent(args.id)
		control.shutdown()

	def experiment(self, args):
		control = ShepherdControl(self.kafka_config)
		control.lattice_experiment(args.number)
		control.shutdown()

	def execute(self):
		args = self.args
		if args.command == 'lattice':
			self.lattice(args)

		elif args.command == 'ecoli':
			self.ecoli(args)

		else:
			super(EnvironmentCommand, self).execute()


if __name__ == '__main__':
	command = EnvironmentCommand()
	command.execute()
