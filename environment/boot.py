from __future__ import absolute_import, division, print_function

import errno
import os
import argparse

from agent.outer import Outer
from agent.inner import Inner

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
	def __init__(self, kafka_config):
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
		self.outer = EnvironmentAgent('EnvironmentSpatialLattice', kafka_config, self.environment)


class BootEcoli(object):
	'''
	This class initializes an EcoliSimulation, passes it to the `Inner` agent, and launches the simulation.
	The EcoliSimulation is initialized by passing it a pathname to sim_data, along with simulation parameters.
	'''
	def __init__(self, agent_id, kafka_config, working_dir,
			variant_type='wildtype', variant_index=0, seed=0):
		self.agent_id = agent_id

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
			self.simulation)


def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	# TODO (eran) share argparse code with agent/boot.py
	# One way to do that is via a base class similar to ScriptBase.py or shared subroutines.
	# another way is argparse.ArgumentParser(parents=[parent_parser]): https://docs.python.org/2/library/argparse.html?highlight=argparse#parents

	parser = argparse.ArgumentParser(
		description='Run an agent for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=['ecoli', 'lattice'],
		help='which command to boot')

	parser.add_argument(
		'--id',
		help='unique identifier for a new simulation agent')

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
		'subscribe_topics': []}

	if args.command == 'lattice':
		BootEnvironmentSpatialLattice(kafka_config)

	elif args.command == 'ecoli':
		if not args.id:
			raise ValueError('the "ecoli" command needs an --id argument')

		BootEcoli(
			args.id, kafka_config, args.working_dir,
			variant_type=args.variant, variant_index=args.index, seed=args.seed
			)

if __name__ == '__main__':
	main()
