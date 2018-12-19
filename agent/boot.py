from __future__ import absolute_import, division, print_function

import json
import argparse

from agent.control import DEFAULT_KAFKA_CONFIG
from agent.outer import Outer
from agent.inner import Inner
from agent.stub import SimulationStub, EnvironmentStub


def boot_outer(agent_id, agent_type, agent_config):
	"""
	Initialize the `EnvironmentStub`, pass it to the `Outer` agent and launch the process.

	This is a demonstration of how to initialize an Outer component. In the place of 
	`EnvironmentStub` you would substitute your own environment class that meets the interface
	defined in `Outer`. 
	"""

	volume = 1
	concentrations = {
		'yellow': 5,
		'green': 11,
		'red': 44,
		'blue': 12}

	environment = EnvironmentStub(volume, concentrations)
	print('outer! {}'.format(agent_config))
	return Outer(agent_id, agent_type, agent_config, environment)

def boot_inner(agent_id, agent_type, agent_config):
	"""
	Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

	This is a demonstration of how to initialize an Inner component. When creating your 
	own simulation you would supply a class that meets the same interface as the `SimulationStub`
	that would be driven by the Inner agent in response to messages from its corresponding 
	Outer agent.
	"""
	if 'outer_id' not in agent_config:
		raise ValueError("--outer-id required")

	agent_id = agent_id
	outer_id = agent_config['outer_id']
	simulation = SimulationStub()

	print('inner... {}'.format(agent_config))
	return Inner(
		agent_id,
		outer_id,
		agent_type,
		agent_config,
		simulation)

class BootAgent(object):
	"""
	Boot agents from the command line.
	"""

	def __init__(self):
		self.description='Boot agents for the environmental context simulation'
		self.agent_types = {
			'outer': boot_outer,
			'inner': boot_inner}

	def add_arguments(self, parser):
		parser.add_argument(
			'--id',
			required=True,
			help='id of the new agent')

		parser.add_argument(
			'--outer-id',
			help="ID of the new agent's outer environment agent")

		parser.add_argument(
			'--type',
			required=True,
			choices=self.agent_types,
			help='type of the new agent')

		parser.add_argument(
			'--config',
			default='{}',
			help='''JSON configuration dictionary for the new agent.''')

		return parser

	def execute(self):
		parser = argparse.ArgumentParser(description=self.description)
		parser = self.add_arguments(parser)
		parse_args = parser.parse_args()

		args = vars(parse_args)
		agent_config = dict(json.loads(parse_args.config))
		agent_config.setdefault('kafka_config', DEFAULT_KAFKA_CONFIG)
		if args['outer_id']:
			agent_config.setdefault('outer_id', args['outer_id'])

		boot = self.agent_types[args['type']]
		agent = boot(args['id'], args['type'], agent_config)
		agent.start()

if __name__ == '__main__':
	boot = BootAgent()
	boot.execute()
