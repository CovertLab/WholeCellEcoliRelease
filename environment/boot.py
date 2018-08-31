from __future__ import absolute_import, division, print_function

import os
import cPickle
import argparse

from agent.outer import Outer
from agent.inner import Inner

from environment.two_dim_lattice import EnvironmentSpatialLattice

# Raw data class
from reconstruction.ecoli.knowledge_base_raw import KnowledgeBaseEcoli

from models.ecoli.sim.simulation import EcoliSimulation

from wholecell.fireworks.firetasks import VariantSimDataTask



class BootEnvironmentSpatialLattice(object):
	def __init__(self, kafka_config):
		raw_data = KnowledgeBaseEcoli()
		# create a dictionary with all saved environments
		self.environment_dict = {}
		for label in dir(raw_data.condition.environment):
			if label.startswith("__"):
				continue
			self.environment_dict[label] = {}
			# initiate all molecules with 0 concentrations
			for row in raw_data.condition.environment_molecules:
				self.environment_dict[label].update({row["molecule id"]: 0}) #* (units.mmol / units.L)})
			# update non-zero concentrations
			molecule_concentrations = getattr(raw_data.condition.environment, label)
			for row in molecule_concentrations:
				self.environment_dict[label].update({row["molecule id"]: row["concentration"].asNumber()})

		# TODO (Eran) don't hardcode initial environment, get this from timeseries
		concentrations = self.environment_dict['minimal']

		self.environment = EnvironmentSpatialLattice(concentrations)
		self.outer = Outer(str(self.environment.agent_id), kafka_config, self.environment)


class BootEcoli(object):
	'''
	This class initializes an EcoliSimulation, passes it to the `Inner` agent, and launches the simulation.
	The EcoliSimulation is initialized by passing it directions to sim_data, along with hardcoded simulation parameters.
	'''
	def __init__(self, agent_id, kafka_config, working_dir):
		self.agent_id = agent_id

		sim_data_fit = '{}/out/manual/kb/simData_Most_Fit.cPickle'.format(working_dir)
		sim_data_variant = '{}/out/manual/kb/simData_Modified.cPickle'.format(working_dir)
		variant_metadata = '{}/out/manual/metadata'.format(working_dir)
		output_dir = '{}/out/manual/sim_{}/simOut'.format(working_dir, self.agent_id)

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


def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=['ecoli', 'lattice'],
		help='which command to boot')

	parser.add_argument(
		'--id',
		help='unique identifier for simulation agent')

	parser.add_argument(
		'--kafka-host',
		default='127.0.0.1:9092',
		help='address for Kafka server')

	parser.add_argument(
		'--environment-control',
		default='environment_control',
		help='topic the environment will receive control messages on')

	parser.add_argument(
		'--simulation-receive',
		default='environment_broadcast',
		help='topic the simulations will receive messages on')

	parser.add_argument(
		'--simulation-send',
		default='environment_listen',
		help='topic the simulations will send messages on')

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
		'subscribe_topics': []}

	if args.command == 'lattice':
		outer = BootEnvironmentSpatialLattice(kafka_config)

	elif args.command == 'ecoli':
		if not args.id:
			raise ValueError('--id must be supplied for inner command')

		inner = BootEcoli(args.id, kafka_config, args.working_dir)

if __name__ == '__main__':
	main()
