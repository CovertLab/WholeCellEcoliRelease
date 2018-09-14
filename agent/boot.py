from __future__ import absolute_import, division, print_function

import os
import argparse

import agent.event as event
from agent.agent import Agent
from agent.outer import Outer
from agent.inner import Inner
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

	def __init__(self, kafka_config):
		volume = 1
		concentrations = {
			'yellow': 5,
			'green': 11,
			'red': 44,
			'blue': 12}

		self.environment = EnvironmentStub(volume, concentrations)
		self.outer = Outer('EnvironmentStub', kafka_config, self.environment)

class BootInner(object):

	"""
	Initialize the `SimulationStub`, pass it to the `Inner` agent and launch the process.

	This is a demonstration of how to initialize an Inner component. When creating your 
	own simulation you would supply a class that meets the same interface as the `SimulationStub`
	that would be driven by the Inner agent in response to messages from its corresponding 
	Outer agent.
	"""

	def __init__(self, agent_id, kafka_config):
		self.agent_id = agent_id
		self.simulation = SimulationStub()
		self.inner = Inner(
			kafka_config,
			self.agent_id,
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

	def trigger_execution(self):
		self.send(self.kafka_config['environment_control'], {
			'event': event.TRIGGER_EXECUTION})

	def pause_execution(self):
		self.send(self.kafka_config['environment_control'], {
			'event': event.PAUSE_ENVIRONMENT})

	def shutdown_environment(self):
		self.send(self.kafka_config['environment_control'], {
			'event': event.SHUTDOWN_ENVIRONMENT})

	def shutdown_simulation(self, agent_id):
		self.send(self.kafka_config['simulation_receive'], {
			'event': event.SHUTDOWN_SIMULATION,
			'inner_id': agent_id})

def main():
	"""
	Parse the arguments for the command line interface to the simulation and launch the
	respective commands.
	"""

	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'command',
		choices=['inner', 'outer', 'trigger', 'pause', 'shutdown'],
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

	if args.command == 'inner':
		if not args.id:
			raise ValueError('--id must be supplied for inner command')

		BootInner(args.id, kafka_config)

	elif args.command == 'outer':
		BootOuter(kafka_config)

	elif args.command == 'trigger':
		control = EnvironmentControl('environment_control', kafka_config)
		control.trigger_execution()
		control.shutdown()

	elif args.command == 'pause':
		control = EnvironmentControl('environment_control', kafka_config)
		control.pause_execution()
		control.shutdown()

	elif args.command == 'shutdown':
		control = EnvironmentControl('environment_control', kafka_config)

		if not args.id:
			control.shutdown_environment()
		else:
			control.shutdown_simulation(args.id)
		control.shutdown()

if __name__ == '__main__':
	main()
