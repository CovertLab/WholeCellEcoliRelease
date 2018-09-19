from __future__ import absolute_import, division, print_function

import os
import uuid
import numpy as np
import argparse

import agent.event as event
from agent.agent import Agent
from agent.outer import Outer
from agent.inner import Inner
from agent.shepherd import AgentShepherd
from agent.stub import SimulationStub, EnvironmentStub


DEFAULT_KAFKA_CONFIG = {
	'host': '127.0.0.1:9092',
	'simulation_send': 'environment_listen',
	'simulation_receive': 'environment_broadcast',
	'environment_control': 'environment_control',
	'subscribe_topics': []}

class BootOuter(object):

	"""
	Initialize the `EnvironmentStub`, pass it to the `Outer` agent and launch the process.

	This is a demonstration of how to initialize an Outer component. In the place of 
	`EnvironmentStub` you would substitute your own environment class that meets the interface
	defined in `Outer`. 
	"""

	def __init__(self, agent_id, agent_config):
		volume = 1
		concentrations = {
			'yellow': 5,
			'green': 11,
			'red': 44,
			'blue': 12}

		self.environment = EnvironmentStub(volume, concentrations)
		self.outer = Outer(agent_id, agent_config['kafka_config'], self.environment)

class BootInner(object):

	"""
	Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

	This is a demonstration of how to initialize an Inner component. When creating your 
	own simulation you would supply a class that meets the same interface as the `SimulationStub`
	that would be driven by the Inner agent in response to messages from its corresponding 
	Outer agent.
	"""

	def __init__(self, agent_id, agent_config):
		self.agent_id = agent_id
		self.outer_id = agent_config['outer_id']
		self.simulation = SimulationStub()
		self.inner = Inner(
			agent_config['kafka_config'],
			self.agent_id,
			self.outer_id,
			self.simulation)


class EnvironmentControl(Agent):

	"""
	Send messages to the other agents in the system to trigger execution and/or shutdown
	the Outer agent (which sends messages to shutdown all the associated Inner agents) or
	shutdown specific Inner agents directly (which then report back to the Outer agent and
	then terminate).
	"""

	def __init__(self, agent_id, kafka_config=None):
		if kafka_config is None:
			kafka_config = DEFAULT_KAFKA_CONFIG.copy()
		super(EnvironmentControl, self).__init__(agent_id, kafka_config)

	def trigger_execution(self, agent_id):
		self.send(self.kafka_config['environment_control'], {
			'event': event.TRIGGER_EXECUTION,
			'agent_id': agent_id})

	def pause_execution(self, agent_id):
		self.send(self.kafka_config['environment_control'], {
			'event': event.PAUSE_ENVIRONMENT,
			'agent_id': agent_id})

	def shutdown_agent(self, agent_id):
		self.send(self.kafka_config['agent_receive'], {
			'event': event.SHUTDOWN_AGENT,
			'agent_id': agent_id})

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

	def add_inner(self, outer_id, agent_config):
		agent_config['outer_id'] = outer_id
		self.add_agent(
			str(uuid.uuid1()),
			'inner',
			agent_config)

	def add_outer(self, agent_id, agent_config):
		self.add_agent(
			agent_id,
			'outer',
			agent_config)

	def stub_experiment(self, inner_number):
		outer_id = str(uuid.uuid1())
		self.add_outer(outer_id, {})
		for index in range(inner_number):
			self.add_inner(outer_id, {})

def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=[
			'inner',
			'outer',
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
		help='unique identifier for simulation agent')

	parser.add_argument(
		'--number',
		type=int,
		default=3,
		help='number of cell agents to spawn in experiment')

	parser.add_argument(
		'--kafka-host',
		default='127.0.0.1:9092',
		help='address for Kafka server')

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
		help='the directory containing the project files'
	)

	args = parser.parse_args()
	kafka_config = {
		'host': args.kafka_host,
		'agent_receive': args.agent_receive,
		'environment_control': args.environment_control,
		'simulation_receive': args.simulation_receive,
		'simulation_send': args.simulation_send,
		'shepherd_control': args.shepherd_control,
		'environment_visualization': args.environment_visualization,
		'subscribe_topics': []}

	if args.command == 'inner':
		if not args.id:
			raise ValueError('--id must be supplied for inner command')

		BootInner(args.id, {'kafka_config': kafka_config})

	elif args.command == 'outer':
		BootOuter(args.id, {'kafka_config': kafka_config, 'outer_id': args.outer_id})

	elif args.command == 'trigger':
		control = EnvironmentControl('environment_control', kafka_config)
		control.trigger_execution(args.id)
		control.shutdown()

	elif args.command == 'pause':
		control = EnvironmentControl('environment_control', kafka_config)
		control.pause_execution(args.id)
		control.shutdown()

	elif args.command == 'shepherd':
		initializers = {}

		def initialize_inner(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = kafka_config
			agent_config['working_dir'] = args.working_dir
			return BootInner(agent_id, agent_config)

		def initialize_outer(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = kafka_config
			return BootOuter(agent_id, agent_config)

		initializers['inner'] = initialize_inner
		initializers['outer'] = initialize_outer

		shepherd = AgentShepherd('agent-shepherd', kafka_config, initializers)

	elif args.command == 'add':
		control = EnvironmentControl('environment_control', kafka_config)
		control.add_inner(args.id, {})
		control.shutdown()

	elif args.command == 'remove':
		control = EnvironmentControl('environment_control', kafka_config)
		control.remove_agent(args.id)
		control.shutdown()

	elif args.command == 'experiment':
		control = EnvironmentControl('environment_control', kafka_config)
		control.stub_experiment(args.number)
		control.shutdown()

	elif args.command == 'shutdown':
		control = EnvironmentControl('environment_control', kafka_config)
		control.shutdown_agent(args.id)
		control.shutdown()

if __name__ == '__main__':
	main()
