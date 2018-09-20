from __future__ import absolute_import, division, print_function

import numpy as np

import agent.event as event
from agent.agent import Agent

class EnvironmentSimulation(object):
	"""Interface for the Outer agent's Environment simulation."""

	def time(self):
		"""Return the current simulation time for the environment."""

	def add_simulation(self, agent_id, state):
		"""Register an inner agent."""

	def remove_simulation(self, agent_id):
		"""Unregister an inner agent."""

	def simulation_parameters(self, agent_id):
		"""Generate any parameters a simulation may need to know about from the environment,
		such as the current time step.
		"""
		return {}

	def update_from_simulations(self, update, now):
		"""Update the environment's state of the inner agent simulations given the
		changes dictionary mapping agent_id to dictionary of molecule counts.
		"""

	def get_molecule_ids(self):
		"""Return the list of molecule IDs."""

	def get_concentrations(self, now):
		"""Return a dictionary of agent_id to concentrations coming from the environment."""

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

	def __init__(self, agent_id, kafka_config, environment):
		"""
		Construct the Agent.

		Args:
			agent_id (str): Unique identifier for this agent.
			kafka_config (dict): Kafka configuration information with the following keys:
				`host`: the Kafka server host address.
				`environment_control`: The topic this agent will use to listen for trigger
					and shutdown messages.
				`simulation_send`: The topic this agent will use to listen for simulation
					updates from inner agents.
			environment (EnvironmentSimulation): The actual simulation which will perform
				the calculations.
		"""
		self.environment = environment
		self.simulations = {}
		self.paused = True
		self.shutting_down = False

		kafka_config['subscribe_topics'] = [
			kafka_config['simulation_send'],
			kafka_config['environment_control']]

		super(Outer, self).__init__(agent_id, kafka_config)

	def initialize(self):
		print('environment started')

	def finalize(self):
		print('environment shutting down')

	def initialize_simulation(self, message):
		inner_id = message['inner_id']
		environment_time = self.environment.time()
		simulation_time = max(environment_time, message.get('time', environment_time))

		self.simulations[inner_id] = {
			'time': simulation_time,
			'message_id': -1,
			'last_message_id': -1,
			'state': message['state']}

		self.send(self.kafka_config['simulation_receive'], {
			'event': event.SYNCHRONIZE_SIMULATION,
			'inner_id': inner_id,
			'outer_id': self.agent_id,
			'state': self.environment.simulation_parameters(inner_id)})

		self.environment.add_simulation(inner_id, message)

	def update_state(self):
		""" Called before each simulation is updated with the current state of the system. """
		pass

	def send_concentrations(self, now, run_until):
		""" Send updated concentrations to each inner agent. """

		concentrations = self.environment.get_concentrations(now)
		self.update_state()

		for inner_id, simulation in self.simulations.iteritems():
			if inner_id in concentrations:
				simulation['message_id'] += 1
				self.send(self.kafka_config['simulation_receive'], {
					'event': event.ENVIRONMENT_UPDATED,
					'outer_id': self.agent_id,
					'inner_id': inner_id,
					'message_id': simulation['message_id'],
					'concentrations': concentrations[inner_id],
					'run_until': run_until})

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

		This will check first to see if all simulations are ready to advance, then it
		will check to see how long each simulation actually ran, only advancing to
		the earliest time point a simulation hit. In this way the environment is always
		running behind the simulations, and each simulation only runs once all other
		simulations have caught up to it (including the environment).
		"""

		if not self.paused and self.ready_to_advance():
			if self.shutting_down:
				self.send_shutdown()
			else:
				# compare the length of each simulation's run
				ran = np.sort([
					simulation['time']
					for simulation in self.simulations.values()])

				# find the earliest time a simulation ran to
				now = ran[0] if ran.size > 0 else 0

				# find any other (longer) run times
				later = ran[ran > now]

				# apply all the updates received from the simulations to the
				# environment's original time point
				self.environment.update_from_simulations(self.simulations, now)

				# run the environment to the current time point
				self.environment.run_incremental(now)

				# find the next time for simulations to achieve
				run_until = self.environment.time() + self.environment.run_for

				# unless there is an earlier time a simulation arrived at
				if later.size > 0:
					run_until = later[0] 

				print('============= environment | ran: {}, now: {}, later: {}, run_until: {}, time: {}'.format(ran, now, later, run_until, self.environment.time()))

				self.send_concentrations(now, run_until)

	def send_shutdown(self):
		for inner_id, simulation in self.simulations.iteritems():
			self.send(self.kafka_config['simulation_receive'], {
				'outer_id': self.agent_id,
				'inner_id': inner_id,
				'event': event.SHUTDOWN_AGENT})

	def receive(self, topic, message):
		"""
		Receive messages from associated inner agents.

		The environment receives messages from both its associated inner agents and also
		the control agent.

		Control messages:

		* TRIGGER_EXECUTION: Send messages to all registered inner agents to begin execution.
		* SHUTDOWN_ENVIRONMENT: Send messages to inner agents notifying them that the outer agent
		    is shutting down, and wait for acknowledgement before exiting.

		Simulation messages:

		* SIMULATION_INITIALIZED: Registers inner agents that will be driven once the
		    TRIGGER_EXECUTION event is received.
		* SIMULATION_ENVIRONMENT: Received from each inner agent when it has computed its 
		    environmental changes up to the specified `run_until`. The outer agent will wait
		    until it has heard from each simulation, integrate their changes and then calculate
		    the new local environment for each inner agent and respond with an `ENVIRONMENT_UPDATED`
		    message.
		* SIMULATION_SHUTDOWN: Received when an inner agent has completed. Once all inner agents
		    have reported back that they have shut down the outer agent can complete.
		"""

		if message.get('outer_id', message.get('agent_id')) == self.agent_id:
			print('--> {}: {}'.format(topic, message))

			if message['event'] == event.SIMULATION_INITIALIZED:
				self.initialize_simulation(message)

			elif message['event'] == event.TRIGGER_EXECUTION:
				self.paused = False
				self.advance()

			elif message['event'] == event.SIMULATION_ENVIRONMENT:
				if message['inner_id'] in self.simulations:
					simulation = self.simulations[message['inner_id']]

					if message['message_id'] == simulation['message_id']:
						simulation['state'] = message['state']
						simulation['time'] = message['time']
						simulation['last_message_id'] = message['message_id']

						self.advance()

			elif message['event'] == event.PAUSE_ENVIRONMENT:
				self.paused = True

			elif message['event'] == event.SHUTDOWN_AGENT:
				if len(self.simulations) > 0:
					if self.ready_to_advance():
						self.send_shutdown()
					else:
						self.shutting_down = True
				else:
					self.shutdown()

			elif message['event'] == event.SIMULATION_SHUTDOWN:
				if message['inner_id'] in self.simulations:
					gone = self.simulations.pop(message['inner_id'], {'inner_id': -1})
					self.environment.remove_simulation(message['inner_id'])

					print('simulation shutdown: ' + str(gone))

					if not self.simulations:
						self.shutdown()

			else:
				print('unexpected event {}: {}'.format(message['event'], message))
