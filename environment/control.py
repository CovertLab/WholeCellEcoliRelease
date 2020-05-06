from __future__ import absolute_import, division, print_function

import copy
import time
import uuid

from vivarium.environment.control import ShepherdControl, EnvironmentCommand
from wholecell.utils import filepath

class EcoliCommand(EnvironmentCommand):
	"""
	Extend `EnvironmentCommand` with new commands related to the lattice and ecoli experiments
	"""

	def __init__(self):
		description = '''
	`ecoli-experiment [--number N] [--type T] [--working-dir D]` ask the Shepherd to run
		a lattice environment with N agents of type T,
	'''
		choices = ['ecoli-experiment']

		super(EcoliCommand, self).__init__(
			choices=choices,
			description=description)

	def ecoli_experiment(self, args):
		self.require(args, 'number', 'working_dir')

		# define experiment: environment type and agent type
		experiment_id = 'lattice_experiment'
		environment_type = 'lattice'
		agents = {
			'ecoli': 1}

		# overwrite default environment config
		lattice_config = {
			'name': 'lattice_experiment',
			'description': (
				'wcecoli in lattice environment')}

		actor_config = copy.deepcopy(self.actor_config)
		actor_config['boot'] = ['python', '-u', '../wcEcoli/environment/boot.py']

		exp_config = {
			'default_experiment_id': experiment_id,
			'lattice_config': lattice_config,
			'environment_type': environment_type,
			'actor_config': actor_config,
			'agents': agents}

		control = ShepherdControl({'kafka_config': self.get_kafka_config()})
		# args['agent_boot'] = ['python', '-u', '../wcEcoli/environment/boot.py']

		control.init_experiment(args, exp_config)
		control.shutdown()

	def add_arguments(self, parser):
		parser = super(EcoliCommand, self).add_arguments(parser)
		return parser

def run():
	command = EcoliCommand()
	command.execute()

if __name__ == '__main__':
	run()
