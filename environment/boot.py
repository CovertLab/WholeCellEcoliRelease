import json
import uuid
import argparse

from environment.agent import Agent
from environment.outer import Outer
from environment.inner import Inner
from environment.stub import SimulationStub

default_kafka_config = {
	'host': 'localhost:9092',
	'simulation_send': 'environment_listen',
	'simulation_receive': 'environment_broadcast',
	'environment_control': 'environment_control',
	'subscribe_topics': []}

class BootOuter(object):
	def __init__(self, kafka):
		self.outer = Outer(
			kafka,
			['yellow', 'green', 'red', 'blue'],
			1,
			{'yellow': 5,
			 'green': 11,
			 'red': 44,
			 'blue': 12})

class BootInner(object):
	def __init__(self, id, kafka):
		self.id = id
		self.simulation = SimulationStub()
		self.inner = Inner(
			self.id,
			self.simulation,
			kafka)

class EnvironmentControl(Agent):
	def __init__(self, kafka=default_kafka_config):
		id = 'environment_control'
		super(EnvironmentControl, self).__init__(id, kafka)

	def trigger_execution(self):
		self.send(self.kafka['environment_control'], {
			'event': 'TRIGGER_EXECUTION'})

	def shutdown_environment(self):
		self.send(self.kafka['environment_control'], {
			'event': 'SHUTDOWN_ENVIRONMENT'})

	def shutdown_simulation(self, id):
		self.send(self.kafka['simulation_receive'], {
			'event': 'SHUTDOWN_SIMULATION',
			'id': id})

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='Boot the various agents for the environmental context simulation')

	parser.add_argument(
		'role',
		help='which role to boot')

	parser.add_argument(
		'--id',
		default='None',
		help='unique identifier for simulation agent')

	parser.add_argument(
		'--kafka-host',
		default='localhost:9092',
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

	args = parser.parse_args()
	kafka = {
		'host': args.kafka_host,
		'environment_control': args.environment_control,
		'simulation_receive': args.simulation_receive,
		'simulation_send': args.simulation_send,
		'subscribe_topics': []}

	if args.role == 'inner':
		inner = BootInner(args.id, kafka)

	elif args.role == 'outer':
		outer = BootOuter(kafka)

	elif args.role == 'trigger':
		control = EnvironmentControl(kafka)
		control.trigger_execution()
		control.shutdown()

	elif args.role == 'shutdown':
		control = EnvironmentControl(kafka)
		if args.id == 'None':
			control.shutdown_environment()
		else:
			control.shutdown_simulation(args.id)
		control.shutdown()
