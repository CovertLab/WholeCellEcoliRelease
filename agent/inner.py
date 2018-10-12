from __future__ import absolute_import, division, print_function

import uuid

import agent.event as event
from agent.agent import Agent


class CellSimulation(object):
	"""Interface for the Inner agent's cell simulation."""

	def time(self):
		"""Return the current time according to this CellSimulation."""

	def synchronize_state(self, state):
		"""Receive any state from the environment, like current time step."""

	def apply_outer_update(self, update):
		"""Apply the update received from the environment to this simulation."""

	def run_incremental(self, run_until):
		"""Run this CellSimulation until the given time."""

	def generate_inner_update(self):
		"""
		Generate the update that will be sent to the environment based on changes calculated
		by the CellSimulation during `run_incremental(run_until)`.

		If the dictionary returned by this function contains a `division` key it will trigger
		preparations for cell division in the environment. The value for this key is a pair of
		dictionaries, of which the only required key is `id` containing the new daughter id.
		Additional keys are specific to the particular implementation of `EnvironmentSimulation`
		that receives this update (see `outer.py`).
		"""

	def divide(self):
		"""Perform cell division on the simulation and return information about the daughters."""
		return []

	def finalize(self):
		"""Release any resources and perform any final cleanup."""


class Inner(Agent):

	"""
	Inner: an independent cell simulation in a larger environmental context.

	This class wraps an instance of CellSimulation into an Agent and mediates
	the message passing communication with the coordinating Outer agent running
	an environmental simulation.
	"""

	def __init__(self, agent_id, outer_id, agent_type, agent_config, simulation):
		"""
		Construct the agent.

		Args:
		    agent_id (str): Unique identifier for this agent.
		        This agent will only respond to messages addressed to its inner agent_id.
		    outer_id (str): Unique identifier for the outer agent this agent belongs to.
		        All messages to an outer agent will be addressed to this id.
		    agent_type (str): The type of this agent, for coordination with the agent shepherd.
		    agent_config (dict): A dictionary containing any information needed to run this
		        outer agent. The only required key is `kafka_config` containing Kafka configuration
		        information with the following keys:

		        * `host`: the Kafka server host address.
		        * `topics`: a dictionary mapping topic roles to specific topics used by the agent
		            to communicate with other agents. The relevant ones to this agent are:

		            * `cell_receive`: The topic this agent will receive messages on from the 
		                environment or relevant control processes.
		            * `environment_receive`: The topic this agent will send messages to its 
		                associated outer agent (given by `outer_id`) and environmental simulation.
		            * `shepherd_receive`: The topic this agent will send messages on for 
		                adding agents to and removing agents from the environment.
		    simulation (CellSimulation): The actual simulation which will perform the
		        calculations.
		"""

		kafka_config = agent_config['kafka_config']
		kafka_config['subscribe'].append(
			kafka_config['topics']['cell_receive'])

		super(Inner, self).__init__(agent_id, agent_type, agent_config)

		self.outer_id = outer_id
		self.simulation = simulation

	def initialize(self):
		"""
		Initialization: Register this inner agent with the outer agent.
		"""

		now = self.simulation.time()
		state = self.simulation.generate_inner_update()

		self.send(self.topics['environment_receive'], {
			'time': now,
			'event': event.CELL_INITIALIZE,
			'outer_id': self.outer_id,
			'inner_id': self.agent_id,
			'agent_config': self.agent_config,
			'state': state})

	def cell_exchange(self, message):
		"""
		Notify the environment simulation about what changes this cell simulation has produced.
		Also, handle the case of the cell dividing during this update by adding the information
		to its environment update and calling `divide_cell`.
		"""

		stop = self.simulation.time()
		update = self.simulation.generate_inner_update()
		division = update.get('division', [])

		for daughter in division:
			if not 'id' in daughter:
				daughter['id'] = str(uuid.uuid1())

		self.send(self.topics['environment_receive'], {
			'event': event.CELL_EXCHANGE,
			'time': stop,
			'outer_id': self.outer_id,
			'inner_id': self.agent_id,
			'message_id': message['message_id'],
			'state': update})

		if division:
			self.divide_cell({}, division)

	def environment_update(self, message):
		"""
		Apply the update from the environment to the cell simulation and respond with changes.
		"""

		self.simulation.apply_outer_update(message['state'])
		self.simulation.run_incremental(message['run_until'])
		self.cell_exchange(message)

	def divide_cell(self, message, division):
		"""
		Perform agent cell division.

		This sends three messages to the agent shepherd: one `ADD_AGENT` for each new daughter cell,
		and finally a `REMOVE_AGENT` for itself. These new agents will initialize and notify the 
		outer agent, inheriting properties from their parent cell.
		"""

		daughter_ids = []

		for daughter in division:
			agent_id = daughter.get('id', str(uuid.uuid1()))
			daughter_ids.append(agent_id)

			agent_type = daughter.get(
				'type',
				message.get(
					'daughter_type',
					self.agent_type))

			agent_config = dict(
				daughter,
				parent_id=self.agent_id,
				outer_id=self.outer_id)

			self.send(self.topics['shepherd_receive'], {
				'event': event.ADD_AGENT,
				'agent_id': agent_id,
				'agent_type': agent_type,
				'agent_config': agent_config})

		self.send(self.topics['shepherd_receive'], {
			'event': event.REMOVE_AGENT,
			'agent_id': self.agent_id})

	def cell_shutdown(self, message):
		self.send(self.topics['environment_receive'], {
			'event': event.CELL_SHUTDOWN,
			'outer_id': self.outer_id,
			'inner_id': self.agent_id})

		self.shutdown()

	def finalize(self):
		"""
		Trigger any clean up the simulation needs to perform before exiting.
		"""

		self.simulation.finalize()

	def receive(self, topic, message):
		"""
		Respond to messages from the environment.

		The inner agent responds to four messages:

		* ENVIRONMENT_SYNCHRONIZE: Receive any relevant information from the environment before
		    the main cycle of mutual updates begins.
		* ENVIRONMENT_UPDATE: Receive the latest state from the environment simulation. The 
		    relevant keys in this update are:

		    * `state`: a dictionary containing the updated state from the environment.
		    * `run_until`: how long to run the cell simulation until before reporting back the new 
		        environmental changes.
		    * `message_id`: the id of the message as provided by the outer agent,
		        to be returned as an acknowledgement that the message was processed along with 
		        the updated environmental changes.

		* DIVIDE_CELL: Perform cell division immediately, whether the simulation is ready or not.
		* SHUTDOWN_AGENT: Shutdown this agent and terminate the process.
		"""

		if message.get('inner_id', message.get('agent_id')) == self.agent_id:
			print('--> {}: {}'.format(topic, message))

			if message['event'] == event.ENVIRONMENT_SYNCHRONIZE:
				self.simulation.synchronize_state(message['state'])

			elif message['event'] == event.ENVIRONMENT_UPDATE:
				self.environment_update(message)

			elif message['event'] == event.DIVIDE_CELL:
				self.simulation.divide()

			elif message['event'] == event.PAUSE_AGENT:
				pass

			elif message['event'] == event.SHUTDOWN_AGENT:
				self.cell_shutdown(message)

			else:
				print('unexpected event {}: {}'.format(message['event'], message))
