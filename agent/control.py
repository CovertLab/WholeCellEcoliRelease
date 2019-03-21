from __future__ import absolute_import, division, print_function

import os
import re
import copy
import json
import uuid
import argparse

import agent.event as event
from agent.agent import Agent

DEFAULT_KAFKA_CONFIG = {
	'host': '127.0.0.1:9092',
	'topics': {
		'agent_receive': 'agent-receive',
		'environment_receive': 'environment-receive',
		'cell_receive': 'cell-receive',
		'shepherd_receive': 'shepherd-receive',
		'visualization_receive': 'environment-state'},
	'subscribe': []}

class AgentControl(Agent):
	"""Send messages to agents in the system to control execution."""

	def __init__(self, agent_id, agent_config=None):
		if 'kafka_config' not in agent_config:
			agent_config['kafka_config'] = copy.deepcopy(DEFAULT_KAFKA_CONFIG)

		super(AgentControl, self).__init__(agent_id, 'control', agent_config)

	def trigger_execution(self, agent_id=''):
		"""Start or resume simulation."""
		if agent_id:
			self.send(self.topics['environment_receive'], {
				'event': event.TRIGGER_AGENT,
				'agent_id': agent_id})
		else:
			self.send(self.topics['shepherd_receive'], {
				'event': event.TRIGGER_ALL})

	def pause_execution(self, agent_id=''):
		if agent_id:
			self.send(self.topics['environment_receive'], {
				'event': event.PAUSE_AGENT,
				'agent_id': agent_id})
		else:
			self.send(self.topics['shepherd_receive'], {
				'event': event.PAUSE_ALL})

	def shutdown_agent(self, agent_id=''):
		if agent_id:
			self.send(self.topics['agent_receive'], {
				'event': event.SHUTDOWN_AGENT,
				'agent_id': agent_id})
		else:
			self.send(self.topics['shepherd_receive'], {
				'event': event.SHUTDOWN_ALL})

	def divide_cell(self, agent_id):
		self.send(self.topics['cell_receive'], {
			'event': event.DIVIDE_CELL,
			'agent_id': agent_id})

	# TODO (Ryan): set this up to send messages to a particular shepherd.
	def add_agent(self, agent_id, agent_type, agent_config):
		kafka = DEFAULT_KAFKA_CONFIG.copy()
		if 'kafka_config' in agent_config:
			kafka.update(agent_config['kafka_config'])
		agent_config['kafka_config'] = kafka

		self.send(self.topics['shepherd_receive'], {
			'event': event.ADD_AGENT,
			'agent_id': agent_id,
			'agent_type': agent_type,
			'agent_config': agent_config})

	def remove_agent(self, agent_query):
		"""
		Remove an agent given either its id or a prefix of its id.

		Args:
		    agent_query (dict): contains either the key `agent_id` or `agent_prefix`.
		        If given an `agent_id`, matches that id exactly.
		        If given `agent_prefix` it will remove any agent whose id has that prefix
                """

		remove = dict(agent_query, event=event.REMOVE_AGENT)
		self.send(self.topics['shepherd_receive'], remove)

	def stub_experiment(self, inner_number):
		outer_id = str(uuid.uuid1())
		self.add_agent(outer_id, 'outer', {})
		for index in range(inner_number):
			self.add_agent(str(uuid.uuid1()), 'inner', {'outer_id': outer_id})

class AgentCommand(object):
	"""
	Control simulations from the command line.

	This class provides a means to send messages to simulations running in a
	distributed environment from the command line. Override `add_arguments` to 
	add more arguments to the argument parser. To respond to more commands,
	supply a `choices` list and implement the methods it names.
	"""

	def __init__(self, choices, description=None):
		self.default_choices = [
			'experiment',
			'add',
			'remove',
			'run',
			'pause',
			'divide',
			'shutdown']
		self.choices = self.default_choices + choices

		if not description:
			description='Boot agents for the environmental context simulation'
		parser = argparse.ArgumentParser(
			description=description,
			formatter_class=argparse.RawDescriptionHelpFormatter)
		self.parser = self.add_arguments(parser)
		self.args = self.parser.parse_args()
		self.kafka_config = {
			'host': self.args.kafka_host,
			'topics': {
				'agent_receive': self.args.agent_receive,
				'environment_receive': self.args.environment_receive,
				'cell_receive': self.args.cell_receive,
				'shepherd_receive': self.args.shepherd_receive,
				'visualization_receive': self.args.visualization_receive},
			'subscribe': []}

	def add_arguments(self, parser):
		parser.add_argument(
			'command',
			choices=self.choices,
			help='which command to run')

		parser.add_argument(
			'--id',
			help='unique identifier for simulation agent')

		parser.add_argument(
			'--type',
			default='ecoli',
			help='type of the agent to add (default: ecoli)')

		parser.add_argument(
			'--config',
			default='{}',
			help='json dictionary of config values to pass to the agent')

		parser.add_argument(
			'--outer-id',
			help='unique identifier for outer agent this inner agent will join')

		parser.add_argument(
			'--prefix',
			help='matches any id with the given prefix')

		parser.add_argument(
			'-n', '--number',
			type=int,
			default=0,
			help='number of cell agents to spawn in the experiment (default: 0)')

		parser.add_argument(
			'--kafka-host',
			default=DEFAULT_KAFKA_CONFIG['host'],
			help='address for Kafka server')

		parser.add_argument(
			'--agent-receive',
			default='agent-receive',
			help='topic agents will receive control messages on')

		parser.add_argument(
			'--environment-receive',
			default='environment-receive',
			help='topic the environment will receive messages on')

		parser.add_argument(
			'--cell-receive',
			default='cell-receive',
			help='topic the simulations will receive messages on')

		parser.add_argument(
			'--shepherd-receive',
			default='shepherd-receive',
			help='topic the shepherd will receive messages on')

		parser.add_argument(
			'--visualization-receive',
			default='environment-state',
			help='topic the environment will send state information on')

		parser.add_argument(
			'--working-dir',
			default=os.getcwd(),
			help='the directory containing the project files')

		return parser

	def require(self, args, *argnames):
		for name in argnames:
			if args.get(name) is None:
				raise ValueError('--{} needed'.format(name))

	def run(self, args):
		control = AgentControl('control', self.kafka_config)
		control.trigger_execution(args['id'])
		control.shutdown()

	def pause(self, args):
		control = AgentControl('control', self.kafka_config)
		control.pause_execution(args['id'])
		control.shutdown()

	def add(self, args):
		self.require(args, 'id', 'type')
		control = AgentControl('control', self.kafka_config)
		config = dict(args['config'], outer_id=args['id'])
		control.add_agent(str(uuid.uuid1()), args['type'] or 'ecoli', config)
		control.shutdown()

	def remove(self, args):
		control = AgentControl('control', self.kafka_config)
		if args['id']:
			control.remove_agent({'agent_id': args['id']})
		elif args['prefix']:
			control.remove_agent({'agent_prefix': args['prefix']})
		else:
			raise ValueError('either --id or --prefix must be provided')
		control.shutdown()

	def divide(self, args):
		self.require(args, 'id')
		control = AgentControl('control', self.kafka_config)
		control.divide_cell(args['id'])
		control.shutdown()

	def experiment(self, args):
		self.require(args, 'number')
		control = AgentControl('control', self.kafka_config)
		control.stub_experiment(args['number'])
		control.shutdown()

	def shutdown(self, args):
		control = AgentControl('control', self.kafka_config)
		control.shutdown_agent(args['id'])
		control.shutdown()

	def execute(self):
		command = self.args.command
		if command in self.choices:
			args = vars(self.args)
			args['config'] = json.loads(self.args.config)
			underscore = re.sub(r'-+', '_', command)
			getattr(self, underscore)(args)
		else:
			raise ValueError("unrecognized command: {}".format(command))

if __name__ == '__main__':
	command = AgentCommand([])
	command.execute()
