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

class AgentCommand(object):
	"""
	Control simulations from the command line.

	This class provides a means to send messages to simulations running in a
	distributed environment from the command line. Override `add_arguments` to 
	add more arguments to the argument parser. Override `execute` to respond to more 
	commands. Override `shepherd_initializers` to provide more types of agents the
	AgentShepherd can spawn.
	"""

	def __init__(self, choices, description=None):
		self.default_choices = [
			'inner',
			'outer',
			'shepherd',
			'experiment',
			'add',
			'remove',
			'trigger',
			'pause',
			'shutdown']
		self.choices = self.default_choices + choices

		if not description:
			description='Boot agents for the environmental context simulation'
		parser = argparse.ArgumentParser(description=description)
		self.parser = self.add_arguments(parser)
		self.args = self.parser.parse_args()
		self.kafka_config = {
			'host': self.args.kafka_host,
			'agent_receive': self.args.agent_receive,
			'environment_control': self.args.environment_control,
			'simulation_receive': self.args.simulation_receive,
			'simulation_send': self.args.simulation_send,
			'shepherd_control': self.args.shepherd_control,
			'environment_visualization': self.args.environment_visualization,
			'subscribe_topics': []}

	def add_arguments(self, parser):
		parser.add_argument(
			'command',
			choices=self.choices,
			help='which command to boot')

		parser.add_argument(
			'--id',
			help='unique identifier for simulation agent')

		parser.add_argument(
			'--outer-id',
			help='unique identifier for outer agent this inner agent will join')

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
			help='the directory containing the project files')

		return parser

	def inner(self, args):
		if not args.id:
			raise ValueError('--id must be supplied for inner command')
		if not args.outer_id:
			raise ValueError('--outer-id must be supplied for inner command')

		BootInner(args.id, {
			'kafka_config': self.kafka_config,
			'outer_id': args.outer_id})

	def outer(self, args):
		if not args.id:
			raise ValueError('--id must be supplied for outer command')

		BootOuter(args.id, {'kafka_config': self.kafka_config})

	def trigger(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.trigger_execution(args.id)
		control.shutdown()

	def pause(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.pause_execution(args.id)
		control.shutdown()

	def shepherd_initializers(self, args):
		initializers = {}

		def initialize_inner(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = self.kafka_config
			agent_config['working_dir'] = args.working_dir
			return BootInner(agent_id, agent_config)

		def initialize_outer(agent_id, agent_config):
			agent_config = dict(agent_config)
			agent_config['kafka_config'] = self.kafka_config
			return BootOuter(agent_id, agent_config)

		initializers['inner'] = initialize_inner
		initializers['outer'] = initialize_outer

		return initializers

	def shepherd(self, args):
		initializers = self.shepherd_initializers(args)
		shepherd = AgentShepherd('agent-shepherd', self.kafka_config, initializers)

	def add(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.add_inner(args.id, {})
		control.shutdown()

	def remove(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.remove_agent(args.id)
		control.shutdown()

	def experiment(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.stub_experiment(args.number)
		control.shutdown()

	def shutdown(self, args):
		control = EnvironmentControl('environment_control', self.kafka_config)
		control.shutdown_agent(args.id)
		control.shutdown()

	def execute(self):
		args = self.args

		if args.command == 'inner':
			self.inner(args)

		elif args.command == 'outer':
			self.outer(args)

		elif args.command == 'trigger':
			self.trigger(args)

		elif args.command == 'pause':
			self.pause(args)

		elif args.command == 'shepherd':
			self.shepherd(args)

		elif args.command == 'add':
			self.add(args)

		elif args.command == 'remove':
			self.remove(args)

		elif args.command == 'experiment':
			self.experiment(args)

		elif args.command == 'shutdown':
			self.shutdown(args)

if __name__ == '__main__':
	command = AgentCommand([])
	command.execute()
