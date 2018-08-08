import json

from environment.agent import Agent

class Inner(Agent):

	"""
	Inner: acts as an independent simulation in a larger environmental context.

	This class wraps a simulation and mediates the communication between messages
	received from the larger environmental context and the operation of the simulation.

	The interface a simulation must provide to work inside of the this class is composed
	of the following methods:

	* simulation.initialize_local_environment()
	    Perfom any setup required for tracking changes to the local environment.

	* simulation.set_local_environment(molecule_ids, concentrations)
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

	If this interface is fulfilled, the simulation can be an Agent in the larger
	environmental simulation. Each Inner Agent will be run in its own thread/process and
	communicate with the larger environmental simulation through message passing.
	"""

	def __init__(self, id, simulation, kafka):
		"""
		Initialize the agent.

		Args:
		    id (string): Unique identifier for this agent in the environmental simulation.
		        When the agent receives messages, it will filter out and respond to only 
		        those containing its `id`.
		    simulation (Simulation): The actual simulation which will perform the calculations.
		    kafka (dict): Kafka configuration information with the following keys:
		        `host`: the Kafka host.
		        `simulation_receive`: The topic the environmental simulation will be sending
		            messages on.
		        `simulation_send`: The topic the environment will be listening to for 
		            updates from the individual simulation agents.
		"""

		self.simulation = simulation
		self.simulation.initialize_local_environment()
		kafka['subscribe_topics'] = [kafka['simulation_receive']]

		super(Inner, self).__init__(id, kafka)

	def initialize(self):
		""" Announce the existence of this simulation agent to the larger environmental context. """

		self.send(self.kafka['simulation_send'], {
			'event': 'SIMULATION_INITIALIZED',
			'id': self.id})

	def finalize(self):
		""" Trigger any clean up the simulation needs to perform before exiting. """

		self.simulation.finalize()

	def receive(self, topic, message):
		"""
		Respond to messages from the environment.

		The Inner Agent responds to only two message: ENVIRONMENT_UPDATED and SHUTDOWN_SIMULATION.
		SHUTDOWN_SIMULATION is called when the system as a whole is shutting down.
		ENVIRONMENT_UPDATED is where the real work of the agent is performed. It receives a 
		message from the environment containing the following keys:

		* `molecule_ids`: ids of the molecules to track.
		* `concentrations`: a dictionary containing the updated local concentrations.
		* `run_until`: how long to run the simulation for before reporting back the new 
		    environmental changes.
		* `message_id`: the id of the message as provided by the environmental simulation,
		    to be returns as an acknowledgement that the message was processed along with 
		    the updated environmental changes.

		Given this, the agent sets the simulations local environment, runs until the given 
		time and responds with a `SIMULATION_ENVIRONMENT` message containing the local changes
		as calculated by the simulation.
		"""

		if message['id'] == self.id:
			print('--> ' + topic + ': ' + str(message))

			if message['event'] == 'ENVIRONMENT_UPDATED':
				self.simulation.set_local_environment(
					message['molecule_ids'],
					message['concentrations'])

				self.simulation.run_incremental(message['run_for'] + self.simulation.time())

				stop = self.simulation.time()
				changes = self.simulation.get_environment_change()

				self.send(self.kafka['simulation_send'], {
					'event': 'SIMULATION_ENVIRONMENT',
					'id': self.id,
					'message_id': message['message_id'],
					'time': stop,
					'changes': changes})

			if message['event'] == 'SHUTDOWN_SIMULATION':
				self.send(self.kafka['simulation_send'], {
					'event': 'SIMULATION_SHUTDOWN',
					'id': self.id})

				self.shutdown()
