from __future__ import absolute_import, division, print_function

import agent.event as event
from agent.agent import Agent


class CellSimulation(object):
	"""Interface for the Inner agent's cell simulation."""

	def time(self):
		"""Return the current time according to this CellSimulation."""

	def initialize_local_environment(self):
		"""Perform any setup required for tracking changes to the local environment."""

	def synchronize_state(self, state):
		"""Receive any state from the environment, like current time step."""

	def set_local_environment(self, concentrations):
		"""Ingest a dictionary of the current chemical concentrations in the
		local environment.
	    """

	def run_incremental(self, run_until):
		"""Run this CellSimulation until the given time."""

	def get_environment_change(self):
		"""Return the accumulated changes to the local environment as calculated
		by the CellSimulation during `run_incremental(run_until)`.
		"""

	def finalize(self):
		"""Release any resources and perform any final cleanup."""


class Inner(Agent):

	"""
	Inner: an independent cell simulation in a larger environmental context.

	This class wraps an instance of CellSimulation into an Agent and mediates
	the message passing communication with the coordinating Outer agent running
	an environmental simulation.
	"""

	def __init__(self, kafka_config, agent_id, outer_id, simulation):
		"""
		Construct the agent.

		Args:
			kafka_config (dict): Kafka configuration information with the following keys:
				`host`: the Kafka server host address.
				`simulation_receive`: The topic this agent will listen to.
				`simulation_send`: The topic this agent will use to send simulation
					updates to the environment.
			agent_id (str): Unique identifier for this agent.
				This agent will only respond to messages addressed to its inner agent_id.
			simulation (CellSimulation): The actual simulation which will perform the
				calculations.
		"""

		self.outer_id = outer_id
		self.simulation = simulation
		self.simulation.initialize_local_environment()
		kafka_config['subscribe_topics'] = [kafka_config['simulation_receive']]

		super(Inner, self).__init__(agent_id, kafka_config)

	def initialize(self):
		"""Initialization: Register this inner agent with the outer agent."""

		now = self.simulation.time()
		state = self.simulation.get_environment_change()

		self.send(self.kafka_config['simulation_send'], {
			'time': now,
			'event': event.SIMULATION_INITIALIZED,
			'outer_id': self.outer_id,
			'inner_id': self.agent_id,
			'state': state})

	def finalize(self):
		""" Trigger any clean up the simulation needs to perform before exiting. """

		self.simulation.finalize()

	def receive(self, topic, message):
		"""
		Respond to messages from the environment.

		The inner agent responds to only two message: ENVIRONMENT_UPDATED and SHUTDOWN_SIMULATION.
		SHUTDOWN_SIMULATION is called when the system as a whole is shutting down.
		ENVIRONMENT_UPDATED is where the real work of the agent is performed. It receives a 
		message from the outer agent containing the following keys:

		* `concentrations`: a dictionary containing the updated local concentrations.
		* `run_until`: how long to run the simulation until before reporting back the new 
		    environmental changes.
		* `message_id`: the id of the message as provided by the outer agent,
		    to be returned as an acknowledgement that the message was processed along with 
		    the updated environmental changes.

		Given this, the agent provides the simulation with the current state of its local
		environment, runs until the given time and responds with a `SIMULATION_ENVIRONMENT`
		message containing the local changes as calculated by the simulation.
		"""

		if message.get('inner_id', message.get('agent_id')) == self.agent_id:
			print('--> {}: {}'.format(topic, message))

			if message['event'] == event.ENVIRONMENT_UPDATED:
				self.simulation.set_local_environment(
					message['concentrations'])

				self.simulation.run_incremental(message['run_until'])

				stop = self.simulation.time()
				update = self.simulation.get_environment_change()

				self.send(self.kafka_config['simulation_send'], {
					'event': event.SIMULATION_ENVIRONMENT,
					'time': stop,
					'outer_id': self.outer_id,
					'inner_id': self.agent_id,
					'message_id': message['message_id'],
					'state': update})

			elif message['event'] == event.SYNCHRONIZE_SIMULATION:
				self.simulation.synchronize_state(message['state'])

			elif message['event'] == event.SHUTDOWN_AGENT:
				self.send(self.kafka_config['simulation_send'], {
					'event': event.SIMULATION_SHUTDOWN,
					'outer_id': self.outer_id,
					'inner_id': self.agent_id})

				self.shutdown()

			else:
				print('unexpected event {}: {}'.format(message['event'], message))

