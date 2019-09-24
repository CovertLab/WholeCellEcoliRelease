from __future__ import absolute_import, division, print_function

import time
import uuid

from lens.environment.control import ShepherdControl, EnvironmentCommand
from wholecell.utils import filepath


class EcoliControl(ShepherdControl):
	"""Send messages to actors in the system to control execution."""

	def __init__(self, agent_config):
		super(EcoliControl, self).__init__(agent_config)


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

		control = EcoliControl({'kafka_config': self.kafka_config})
		args['agent_boot'] = ['python', '-u', '../wcEcoli/environment/boot.py']

		control.lattice_experiment(args)
		control.shutdown()

	def add_arguments(self, parser):
		parser = super(EcoliCommand, self).add_arguments(parser)
		return parser

def run():
	command = EcoliCommand()
	command.execute()

if __name__ == '__main__':
	run()
