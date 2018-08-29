from __future__ import absolute_import, division, print_function

import agent.event as event
from agent.agent import Agent


class Inner(Agent):

	"""
	Inner: acts as an independent simulation in a larger environmental context.

	This class wraps a simulation and mediates the communication between messages
	received from the coordinating outer agent and the operation of the simulation.

	The interface a simulation must provide to work inside of this class is composed
	of the following methods:

	* simulation.initialize_local_environment()
	    Perform any setup required for tracking changes to the local environment.

	* simulation.set_local_environment(concentrations)
	    Receive a set of molecule ids to track and a dictionary containing the current
	    current concentrations in the local environment.

	* simulation.time()
	    Return the current time according to the simulation.

	* simulation.run_incremental(run_until)
	    Run the simulation until the given time.

	* simulation.get_environment_change()
	    Return the accumulated changes to the local environment as calculated by the 
	    simulation during `run_incremental(run_until)`.

	* simulation.finalize()
	    Release any resources and perform any final cleanup.

	If this interface is fulfilled, the simulation can be an agent in the larger
	environmental simulation. Each inner agent will be run in its own thread/process and
	communicate with the outer agent through message passing.
	"""

	def __init__(self, kafka_config, agent_id, simulation):
		"""
		Initialize the agent.

		Args:
		    kafka_config (dict): Kafka configuration information with the following keys:
		        `host`: the Kafka host.
		        `simulation_receive`: The topic the outer agent will be sending messages on.
		        `simulation_send`: The topic the outer agent will be listening to for 
		            updates from the inner agents.
		    agent_id (string): Unique identifier for this agent.
		        When the agent receives messages, it will filter out and respond to only 
		        those containing its `id`.
		    simulation (Simulation): The actual simulation which will perform the calculations.
		"""

		self.simulation = simulation
		self.simulation.initialize_local_environment()
		kafka_config['subscribe_topics'] = [kafka_config['simulation_receive']]

		super(Inner, self).__init__(agent_id, kafka_config)

	def initialize(self):
		""" Announce the existence of this inner agent to the outer agent. """

		self.send(self.kafka_config['simulation_send'], {
			'event': event.SIMULATION_INITIALIZED,
			'inner_id': self.id})

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

		if message['inner_id'] == self.id:
			print('--> {}: {}'.format(topic, message))

			if message['event'] == event.ENVIRONMENT_UPDATED:
				self.simulation.set_local_environment(
					message['concentrations'])

				self.simulation.run_incremental(message['run_until'])

				stop = self.simulation.time()
				changes = self.simulation.get_environment_change()

				self.send(self.kafka_config['simulation_send'], {
					'event': event.SIMULATION_ENVIRONMENT,
					'inner_id': self.id,
					'message_id': message['message_id'],
					'time': stop,
					'changes': changes})

			elif message['event'] == event.SHUTDOWN_SIMULATION:
				self.send(self.kafka_config['simulation_send'], {
					'event': event.SIMULATION_SHUTDOWN,
					'inner_id': self.id})

				self.shutdown()

			else:
				print('unexpected event {}: {}'.format(message['event'], message))
