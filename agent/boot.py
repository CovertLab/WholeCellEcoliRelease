from __future__ import absolute_import, division, print_function

import uuid
import argparse

from agent.agent import Agent
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
	return Outer(agent_id, agent_type, agent_config, environment)

def boot_inner(agent_id, agent_type, agent_config):
	"""
	Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

	This is a demonstration of how to initialize an Inner component. When creating your 
	own simulation you would supply a class that meets the same interface as the `SimulationStub`
	that would be driven by the Inner agent in response to messages from its corresponding 
	Outer agent.
	"""

	agent_id = agent_id
	outer_id = agent_config['outer_id']
	simulation = SimulationStub()

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
		description='Boot agents for the environmental context simulation'
		parser = argparse.ArgumentParser(description=description)
		self.parser = self.add_arguments(parser)
		self.args = self.parser.parse_args()
		self.agent_types = {
			'outer': boot_outer,
			'inner': boot_inner}

	def add_arguments(self, parser):
		parser.add_argument(
			'id',
			default=uuid.uuid1(),
			help='id of the new agent')

		parser.add_argument(
			'type',
			default='inner',
			help='type of the new agent')

		parser.add_argument(
			'config',
			default='{}',
			help='configuration for the new agent')

	def execute(self):
		args = vars(self.args)
		args['config'] = json.loads(self.args.config)
		boot = self.agent_types[args['type']]
		agent = boot(args['id'], args['type'], args['config'])
		agent.start()

if __name__ == '__main__':
	boot = BootAgent()
	boot.execute()
