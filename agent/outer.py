from __future__ import absolute_import, division, print_function

import agent.event as event
from agent.agent import Agent

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

	The interaction with the environmental simulation is mediated through an interface defined
	by the following functions:

	* environment.time()
	    Return the current simulation time for the environment.

	* environment.add_simulation(id)

	* environment.remove_simulation(id)

	* environment.update_from_simulations(changes)
        `changes` is a dictionary of simulation ids to counts

	* environment.run_simulations_until()
	    Returns a dictionary of simulation ids to time points for each simulation to run until.

	* environment.get_molecule_ids()

	* environment.get_concentrations()
	    Returns a dictionary of simulation ids to concentrations coming from the environment.
	"""

	def __init__(self, agent_id, kafka_config, environment):
		self.environment = environment
		self.simulations = {}
		self.shutting_down = False

		kafka_config['subscribe_topics'] = [
			kafka_config['simulation_send'],
			kafka_config['environment_control']]

		super(Outer, self).__init__(agent_id, kafka_config)

	def initialize(self):
		print('environment started')

	def finalize(self):
		print('environment shutting down')

	def initialize_simulation(self, agent_id):
		changes = {molecule: 0 for molecule in self.environment.get_molecule_ids()}
		self.simulations[agent_id] = {
			'time': 0,
			'message_id': -1,
			'last_message_id': -1,
			'changes': changes}

		self.environment.add_simulation(agent_id)

	def send_concentrations(self):
		""" Send updated concentrations to each inner agent. """

		concentrations = self.environment.get_concentrations()
		run_until = self.environment.run_simulations_until()

		for agent_id, simulation in self.simulations.iteritems():
			simulation['message_id'] += 1
			self.send(self.kafka_config['simulation_receive'], {
				'inner_id': agent_id,
				'message_id': simulation['message_id'],
				'event': event.ENVIRONMENT_UPDATED,
				'concentrations': concentrations[agent_id],
				'run_until': run_until[agent_id]})

		self.environment.run_incremental(run_until[self.environment.agent_id])

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

	def simulation_changes(self):
		return {
			agent_id: simulation['changes']
			for agent_id, simulation in self.simulations.iteritems()
			}

	def send_shutdown(self):
		for agent_id, simulation in self.simulations.iteritems():
			self.send(self.kafka_config['simulation_receive'], {
				'inner_id': agent_id,
				'event': event.SHUTDOWN_SIMULATION})

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

		print('--> {}: {}'.format(topic, message))

		if message['event'] == event.SIMULATION_INITIALIZED:
			self.initialize_simulation(message['inner_id'])

		elif message['event'] == event.TRIGGER_EXECUTION:
			self.send_concentrations()

		elif message['event'] == event.SIMULATION_ENVIRONMENT:
			if message['inner_id'] in self.simulations:
				simulation = self.simulations[message['inner_id']]

				if message['message_id'] == simulation['message_id']:
					simulation['changes'] = message['changes']
					simulation['time'] = message['time']
					simulation['last_message_id'] = message['message_id']

					if self.ready_to_advance():
						if self.shutting_down:
							self.send_shutdown()
						else:
							changes = self.simulation_changes()
							self.environment.update_from_simulations(changes)
							self.send_concentrations()

		elif message['event'] == event.SHUTDOWN_ENVIRONMENT:
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
