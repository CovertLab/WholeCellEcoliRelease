from __future__ import absolute_import, division, print_function

from collections import defaultdict
import math
import uuid
import numpy as np
import os
import time

import agent.event as event
from agent.agent import Agent
import wholecell.utils.filepath as fp


class EnvironmentSimulation(object):
	"""Interface for the Outer agent's Environment simulation."""

	def time(self):
		"""Return the current simulation time for the environment."""
		return 0

	def add_simulation(self, agent_id, simulation):
		"""Register an inner agent."""

	def remove_simulation(self, agent_id):
		"""Unregister an inner agent."""

	def simulation_parameters(self, agent_id):
		"""
		During initialization of a new simulation, generate any parameters a simulation may
		need to know about from the environment, such as the current time step.
		"""
		return {}

	def simulation_state(self, agent_id):
		"""
		Return the state the environment is tracking about the simulation given by `agent_id`.
		"""
		return {}

	def apply_inner_update(self, update, now):
		"""
		Update the environment's state from only the inner agent simulations that have reached
		but not passed the `now` time point, given the `update` dictionary mapping agent_id to
		the changes calculated during their run.
		"""

	def generate_outer_update(self, now):
		"""
		Return a dictionary of agent_id to updates coming from the environment for
		each agent that has run to `now` but not past.
		"""
		return {}

	def apply_parent_state(self, agent_id, simulation):
		"""
		After cell division, this function is called when a new daughter cell is initialized
		by the environment in order to apply any state the environment was tracking about the 
		parent cell (like location and orientation etc). `simulation` is a dict describing
		the new daughter but here, some fields are just copied from the parent.
		"""

	def run_for_time(self):
		"""Return the length of time simulations should run for this time period."""
		return 0

	def max_time(self):
		"""Return the maximum time for the simulations."""
		return float('inf')

	def run_incremental(self, time):
		"""Run the environment's own simulation until the given time."""


class Outer(Agent):

	"""
	Outer: coordinate the communication between inner agents running their individual simulations
	and the containing environmental context.

	This class represents the larger environmental context for each of the individual simulation
	agents. The general flow is that each inner agent will declare its existence to the 
	outer agent until the outer agent receives a signal from the control process to begin execution.
	Once this signal is received the outer agent will broadcast the local environmental
	concentrations to each inner agent, at which point they will perform their local calculations
	and report back with their updated local environment. Once the outer agent receives
	a message from each inner agent these changes are integrated, and the newly updated environmental
	concentrations will be sent back. This loop will continue until the outer agent receives a 
	message to shutdown, when it will send a message to each inner agent to shutdown and wait for
	acknowledgements, at which point it will shutdown itself.

	Inner agents may also be added and removed while the execution is running without interruption.

	The context environmental simulation is an instance of EnvironmentSimulation.
	"""

	def __init__(self, agent_id, agent_type, agent_config, environment):
		"""
		Construct the Agent.

		Args:
		    agent_id (str): Unique identifier for this agent.
		    agent_type (str): The type of this agent, for coordination with the agent shepherd.
		    agent_config (dict): A dictionary containing any information needed to run this
		        outer agent. The only required key is `kafka_config` containing Kafka configuration
		        information with the following keys:

		        * `host`: the Kafka server host address.
		        * `topics`: a dictionary mapping topic roles to specific topics used by the agent
		            to communicate with other agents. The relevant ones to this agent are:

		            * `environment_receive`: The topic this agent will use to listen for
		                any messages addressed to it, either from related simulations or 
		                from control messages from other processes (such as trigger and shutdown).
		            * `cell_receive`: The topic this agent will use to send messages to its
		                associated cell simulations.
		    environment (EnvironmentSimulation): The actual simulation which will perform
		        the calculations.
		"""

		kafka_config = agent_config['kafka_config']
		kafka_config['subscribe'].append(
			kafka_config['topics']['environment_receive'])

		super(Outer, self).__init__(agent_id, agent_type, agent_config)

		self.environment = environment
		self.simulations = {}
		self.paused = True
		self.shutting_down = False

		# Log the child -> parent relationships for lineage analysis.
		working_dir = agent_config.get('working_dir', os.getcwd())
		output_dir = fp.makedirs(working_dir, 'out', 'manual', agent_id)
		self.lineage_filename = os.path.join(output_dir, 'cell_lineage.json')
		self.lineage = {}

		self.update_state()

	def preinitialize(self):
		print('environment started')

	def finalize(self):
		print('environment shutting down')

	def synchronize_new_cell(self, message):
		'''
		Synchronize clocks with a newly declared/initialized cell and register state info about that cell.
		This gets called when a new cell agent optionally "declares" that it's coming into existence (so
		Lens can add a stand-in cell without waiting for its simulation object to initialize) and when the
		new cell agent reports that it's "initialized" as fully ready.
		'''

		inner_id = message['inner_id']

		simulation_time = self.environment.time()
		if self.simulations:
			latest = max([
				simulation['time']
				for agent_id, simulation
				in self.simulations.iteritems()])
			simulation_time = max(simulation_time, latest)

		simulation = self.simulations.setdefault(inner_id, {})

		simulation.update({
			'time': simulation_time,
			'message_id': -1,
			'last_message_id': -2,
			'state': message['state'],
			'agent_config': message['agent_config']})

		self.environment.add_simulation(inner_id, simulation)

		# lineage tracing
		parent_id = simulation.get('parent_id', '')
		if parent_id:
			self.environment.apply_parent_state(inner_id, simulation)

		if inner_id not in self.lineage:
			self.lineage[inner_id] = parent_id
			fp.write_json_file(self.lineage_filename, self.lineage, indent=2)

		self.update_state()

	def cell_declare(self, message):
		'''
		Synchronize the inner agent with any pertinent state, and trigger its initialization.

		After receiving a CELL_DECLARE message, this function sends an ENVIRONMENT_SYNCHRONIZE
		message to the inner agent to trigger its initialization. The ENVIRONMENT_SYNCHRONIZE
		message includes 'state', which can include variable needed for the cell's initialization.
		'''

		self.synchronize_new_cell(message)
		inner_id = message['inner_id']

		# synchronize state of the new cell
		parameters = self.environment.simulation_parameters(inner_id)
		self.send(self.topics['cell_receive'], {
			'event': event.ENVIRONMENT_SYNCHRONIZE,
			'inner_id': inner_id,
			'outer_id': self.agent_id,
			'state': parameters})

	def cell_initialize(self, message):
		"""
		Prepares the initialization of a new cell simulation.

		First, update the cell  so that subclasses can perform whatever operation they need to
		perform when the state of the simulation changes (such as sending state notifications to
		listening visualizations). Finally, the outer agent is advanced if all associated inner
		agents are ready to advance.
		"""

		self.synchronize_new_cell(message)

		inner_id = message['inner_id']
		simulation = self.simulations[inner_id]

		# print('=== initializing simulation {}'.format(simulation))

		simulation.update({
			'message_id': -1,
			'last_message_id': -1})

		self.advance()

	def update_state(self):
		"""
		Called when the overall state of the environment simulation and its associated cell
		simulations is changed so that external listening processes can be notified.
		"""
		pass

	def send_updates(self, now, run_until):
		"""
		Send updates from the environment simulation out to the inner agents and their associated
		cell simulations and tell them when to `run_until`, if they have stopped at `now`. 
		"""

		update = self.environment.generate_outer_update(now)
		self.update_state()

		for inner_id, simulation in self.simulations.iteritems():
			if inner_id in update:
				simulation['message_id'] += 1
				self.send(self.topics['cell_receive'], {
					'event': event.ENVIRONMENT_UPDATE,
					'outer_id': self.agent_id,
					'inner_id': inner_id,
					'message_id': simulation['message_id'],
					'state': update[inner_id],
					'run_until': run_until})

	def cell_exchange(self, message):
		"""
		Handle messages from inner agents about the behavior of their cell simulations.

		This updates the state for each cell simulation and also handles the case of cell division,
		where the lattice state (location, volume) of the parent cell is copied to daughter cells.

		Also, this (along with `cell_initialize`) is the main juncture where the outer simulation
		is advanced if all cell simulations are ready.
		"""

		agent_id = message['inner_id']
		if agent_id in self.simulations:
			simulation = self.simulations[agent_id]

			# check to make sure the message we just received is in reply to the most recent
			# message we sent to the inner agent
			if message['message_id'] == simulation['message_id']:
				state = message['state']
				simulation['state'] = state
				simulation['time'] = message['time']

				# if the update from the inner agent contains a non-empty `division` list,
				# prepare the state of each impending daughter cell. The `division` list
				# contains info dictionaries for the `EnvironmentSimulation`.
				if state.get('division'):
					parent = self.environment.simulation_state(agent_id)
					for index, daughter in enumerate(state['division']):
						daughter_id = daughter.get('id', str(uuid.uuid1()))
						self.simulations[daughter_id] = dict(
							parent,
							time=simulation['time'],
							parent_id=agent_id,
							index=index,
							message_id=0,
							last_message_id=-1)
						# print('=== daughter: {}'.format(self.simulations[daughter_id]))
				else:
					simulation['last_message_id'] = message['message_id']

				self.advance()

	def ready_to_advance(self):
		"""
		Predicate to determine if the outer agent has heard back from all known inner agents,
		in which case the outer agent can proceed to the next step.
		"""

		# TODO(Ryan): replace this with something that isn't O(n^2)
		ready = True
		for agent_id, simulation in self.simulations.iteritems():
			if simulation['message_id'] > simulation['last_message_id']:
				ready = False
				break

		return ready

	def advance(self):
		"""
		Advance the environment once it has heard back from all registered simulations,
		then send out the newly calculated concentrations to each simulation.

		This will check first to see if all cells are ready to advance, then it
		will check to see how long each cell actually ran, only advancing to
		the earliest time point a cell hit. In this way the environment is always
		running behind the cells, and each cell only runs once all other
		cells and the environment have caught up to it.
		"""

		if not self.paused and self.ready_to_advance():
			if self.shutting_down:
				self.send_shutdown()
			elif self.environment.time() > self.environment.max_time():
				print('Passed max_time {}'.format(self.environment.max_time()))
				self.send_shutdown()
			else:
				# compare the length of each simulation's run
				ran = np.sort([
					math.ceil(simulation['time'])
					for simulation
					in self.simulations.values()])

				# find the earliest time a simulation ran to
				now = ran[0] if ran.size > 0 else 0

				# find any other (longer) run times
				later = ran[ran > now]

				# apply all the updates received from the simulations to the
				# environment's original time point
				self.environment.apply_inner_update(self.simulations, now)

				# run the environment to the current time point
				self.environment.run_incremental(now)

				# find the next time for simulations to achieve
				run_until = self.environment.time() + self.environment.run_for_time()

				# unless there is an earlier time a simulation arrived at
				if later.size > 0:
					run_until = later[0] 

				# print('=== environment | ran: {}, now: {}, later: {}, run_until: {}, time: {}'.format(ran, now, later, run_until, self.environment.time()))

				self.send_updates(now, run_until)

	def cell_shutdown(self, message):
		if message['inner_id'] in self.simulations:
			gone = self.simulations.pop(message['inner_id'], {'inner_id': -1})
			self.environment.remove_simulation(message['inner_id'])
			self.update_state()

			print('simulation shutdown: ' + str(gone))

			if not self.simulations:
				self.shutdown()

	def send_shutdown(self):
		for inner_id, simulation in self.simulations.iteritems():
			self.send(self.topics['cell_receive'], {
				'outer_id': self.agent_id,
				'inner_id': inner_id,
				'event': event.SHUTDOWN_AGENT})

	def shutdown_inner(self, message):
		if len(self.simulations) > 0:
			if self.ready_to_advance():
				self.send_shutdown()
			else:
				self.shutting_down = True
		else:
			self.shutdown()

	def receive(self, topic, message):
		"""
		Receive messages from associated inner agents.

		The environment receives messages from both its associated inner agents and also
		the control agent.

		Control messages:

		* TRIGGER_AGENT: Send messages to all registered inner agents to begin execution.
		* PAUSE_AGENT: Stop sending messages until another TRIGGER_AGENT is received.
		* SHUTDOWN_AGENT: Send messages to inner agents notifying them that the outer agent
		    is shutting down, and wait for acknowledgement before exiting.

		Simulation messages:

		* CELL_DECLARE: Declare an inner agent. [How does this differ from CELL_INITIALIZE?]
		* CELL_INITIALIZE: Registers inner agents that will be driven once the
		    TRIGGER_AGENT event is received.
		* CELL_EXCHANGE: Received from each inner agent when it has computed its 
		    environmental changes up to the specified `run_until`. The outer agent will wait
		    until it has heard from each simulation, integrate their changes and then calculate
		    the new local environment for each inner agent and respond with an `ENVIRONMENT_UPDATE`
		    message.
		* CELL_SHUTDOWN: Received when an inner agent has completed. Once all inner agents
		    have reported back that they have shut down the outer agent can complete.
		"""

		if message.get('outer_id', message.get('agent_id')) == self.agent_id:
			self.print_message(topic, message)

			if message['event'] == event.TRIGGER_AGENT:
				self.paused = False
				self.advance()

			elif message['event'] == event.PAUSE_AGENT:
				self.paused = True
				self.update_state()

			elif message['event'] == event.SHUTDOWN_AGENT:
				self.shutdown_inner(message)

			elif message['event'] == event.CELL_DECLARE:
				self.cell_declare(message)

			elif message['event'] == event.CELL_INITIALIZE:
				self.cell_initialize(message)

			elif message['event'] == event.CELL_EXCHANGE:
				self.cell_exchange(message)

			elif message['event'] == event.CELL_SHUTDOWN:
				self.cell_shutdown(message)

			else:
				print('unexpected event {}: {}'.format(message['event'], message))
