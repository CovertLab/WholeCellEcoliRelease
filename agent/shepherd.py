from __future__ import absolute_import, division, print_function

import multiprocessing as mp
import agent.event as event
from agent.agent import Agent

class AgentShepherd(Agent):

	"""
	AgentShepherd is an agent that spawns, tracks and removes other agent processes in
	a multiprocessing environment. The type of agents it is able to spawn is mediated by
	a dictionary of "initializers" that is passed to it on init. Each key is the name of an agent
	type and each value is a function that takes two arguments, an `agent_id` and an `agent_config`
	dictionary. Each time an agent is added, it is done so by calling one of these initializers.

	AgentShepherd responds to two agent lifecycle messages, ADD_AGENT and REMOVE_AGENT, as well as
	three control messages, TRIGGER_ALL, PAUSE_ALL, and SHUTDOWN_ALL, which apply these messages to
	every agent in its process pool.

	* ADD_AGENT: takes an `agent_id`, an `agent_type` which is the key to an initializer, and
	    `agent_config`, which is passed to the initializer when spawning the new agent. Each
	    agent is a separate process using the python `multiprocessing` library.

	* REMOVE_AGENT: takes an `agent_id` to remove or an `agent_prefix`, which will remove all
	    agents whose ids start with the given string. This enables you to remove several agents
	    at once, or to just supply the beginning of a long agent id string (like a uuid), similar
	    to git's behavior when naming commit hashes.
	"""

	def __init__(self, agent_id, agent_config, agent_initializers):
		"""
		Initialize the AgentShepherd with its id, kafka config and a dictionary of initializers,
		which determine what kind of agents the shepherd is able to spawn.

		Args:
		    agent_id (str): A unique identifier for the new agent.
		    agent_config (dict): A dictionary containing any information needed to run this
		        outer agent. The only required key is `kafka_config` containing Kafka configuration
		        information with the following keys:

		        * `host`: the Kafka server host address.
		        * `topics`: a dictionary mapping topic roles to specific topics used by the agent
		            to communicate with other agents. The relevant ones to this agent are:

		            * `shepherd_receive`: The topic this agent will be listening to for directives
		                to add or remove agents from the simulation.
		            * `agent_receive`: The topic this agent will send control messages on for 
		                running, pausing, or shutting down agents.
		    agent_initializers (dict): This is the set of agents this shepherd will be able to
		        spawn. The values are callables that take two arguments, `agent_id` (str) and
		        `agent_config` (dict), and will be passed to `multiprocessing.Process` to spawn
		        a new agent process.
		"""

		kafka_config = agent_config['kafka_config']
		kafka_config['subscribe'].append(
			kafka_config['topics']['shepherd_receive'])

		super(AgentShepherd, self).__init__(agent_id, 'shepherd', agent_config)

		self.agents = {}
		self.agent_initializers = agent_initializers

	def initialize(self):
		print('agent shepherd waiting')


	# TODO(Ryan): add command to list agent state
	def add_agent(self, agent_id, agent_type, agent_config):
		"""
		Add a new agent for the shepherd to track.

		Args:
		    agent_id (str): Unique identifier for the new agent.
		    agent_type (str): Specifies which type of agent to create.
		        Must be one of the keys of this shepherd's `initializers` dictionary.
		    agent_config (dict): Any parameters the agent needs for initialization. This
		        dictionary will be passed to the initializer on invocation along with the agent_id.
		"""

		initializer = self.agent_initializers.get(agent_type, None)
		if initializer:
			process = mp.Process(
				target=initializer,
				name=agent_id,
				args=(agent_id, agent_type, agent_config))

			self.agents[agent_id] = {
				'process': process,
				'id': agent_id,
				'type': agent_type,
				'config': agent_config}

			process.start()

		else:
			print('agent initializer not found for {}'.format(agent_type))

	def remove_prefix(self, agent_prefix):
		"""
		Remove all agents from the pool with an id containing the given prefix.

		Args:
		    agent_prefix (str): This prefix will match all agent ids that begin with the
		        given string. This way you don't need to type a whole uuid, and can even remove
		        multiple agents at once if their ids are logically grouped by prefix.
		"""

		removing = filter(lambda key: key.startswith(agent_prefix), self.agents.iterkeys())
		print('removing agents {}'.format(removing))

		removed = {}
		for agent_id in removing:
			removed[agent_id] = self.remove_agent(agent_id)

		print('removal complete {}'.format(removing))

		return removed

	def remove_agent(self, agent_id):
		"""
		Remove the agent with the given `agent_id`.
		"""

		if agent_id in self.agents:
			agent = self.agents.pop(agent_id)
			if agent['process'].is_alive():
				self.send(self.topics['agent_receive'], {
					'event': event.SHUTDOWN_AGENT,
					'agent_id': agent_id})

			return agent

	def filter_type(self, agents, agent_type):
		"""
		Find agents only of the given `agent_type`
		"""

		if agent_type:
			agents = filter(
				lambda agent: agent['type'] == agent_type,
				agents)
		return agents

	def agent_control(self, agent_event, agents, agent_type=None):
		"""
		Send a control message of the given `agent_event` type to each agent of the given
		`agent_type`, if not `None`, otherwise to every agent.
		"""

		matching = self.filter_type(agents, agent_type)
		for agent in matching:
			self.send(self.topics['agent_receive'], {
				'event': agent_event,
				'agent_id': agent['id']})

	def receive(self, topic, message):
		"""
		This agent receives two agent lifecycle events: ADD_AGENT and REMOVE_AGENT, who's message
		contents match the arguments to the functions above.

		In addition it receives three control messages: TRIGGER_ALL, PAUSE_ALL, and SHUTDOWN_ALL,
		which apply the corresponding control messages to every agent in their agent pool, or just 
		those of a given type if `agent_type` is providing as well.
		"""

		# TODO (Ryan): Make the shepherd respond to messages addressed to it, rather than every
		#     message of the given type.
		# if message['agent_id'] == self.agent_id:

		self.print_message(topic, message)

		if message['event'] == event.ADD_AGENT:
			self.add_agent(
				message['agent_id'],
				message['agent_type'],
				message['agent_config'])

		elif message['event'] == event.REMOVE_AGENT:
			if 'agent_prefix' in message:
				self.remove_prefix(message['agent_prefix'])
			else:
				self.remove_agent(message['agent_id'])

		elif message['event'] == event.TRIGGER_ALL:
			self.agent_control(event.TRIGGER_AGENT, self.agents.values(), message.get('agent_type'))

		elif message['event'] == event.PAUSE_ALL:
			self.agent_control(event.PAUSE_AGENT, self.agents.values(), message.get('agent_type'))

		elif message['event'] == event.SHUTDOWN_ALL:
			self.agent_control(event.SHUTDOWN_AGENT, self.agents.values(), message.get('agent_type'))
