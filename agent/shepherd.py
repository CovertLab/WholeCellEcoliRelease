from __future__ import absolute_import, division, print_function

import multiprocessing as mp
import agent.event as event
from agent.agent import Agent

class AgentShepherd(Agent):

	"""
	AgentShepherd is an agent that spawns, tracks and removes other agent processes in
	a multiprocessing environment. The type of agents it is able to spawn is mediated by
	a dictionary of "initializers" that is passed to it on init. Each key is the name of an agent
	type and each value is a function that takes two arguments, and `agent_id` and an `agent_config`
	dictionary. Each time an agent is added, it is done so by calling one of these initializers.

	AgentShepherd responds to two messages, ADD_AGENT and REMOVE_AGENT.

	* ADD_AGENT: takes an `agent_id`, an `agent_type` which is the key to an initializer, and
	    `agent_config`, which is passed to the initializer when spawning the new agent. Each
	    agent is a separate process using the python `multiprocessing` library.
	* REMOVE_AGENT: takes an `agent_prefix` which will remove all agents whose ids start with
	    the given string. This enables you to remove several agents at once, or to just supply
	    the beginning of a long agent id string (like a uuid), similar to git's behavior when
	    naming commit hashes.
	"""

	def __init__(self, agent_id, kafka_config, agent_initializers):
		"""
		Initialize the AgentShepherd with its id, kafka config and a dictionary of initializers,
		which determine what kind of agents the shepherd is able to spawn.
		"""

		self.agents = {}
		self.agent_initializers = agent_initializers

		kafka_config['subscribe_topics'] = [kafka_config['shepherd_control']]
		super(AgentShepherd, self).__init__(agent_id, kafka_config)

	def initialize(self):
		print('agent shepherd waiting')

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
			process = mp.Process(target=initializer, args=(agent_id, agent_config))
			process.start()
			self.agents[agent_id] = {
				'process': process,
				'type': agent_type,
				'config': agent_config}
		else:
			print('agent initializer not found for {}'.format(agent_type))

	def remove_agent(self, agent_prefix):
		"""
		Remove an agent from the pool given a prefix of its id.

		Args:
		    agent_prefix (str): This prefix will match all agent ids that begin with the
		        given string. This way you don't need to type a whole uuid, and can even remove
		        multiple agents at once if their ids are logically grouped by prefix.
		"""

		removing = filter(lambda key: key.startswith(agent_prefix), self.agents.iterkeys())
		print('removing agents {}'.format(removing))

		removed = {}
		for key in removing:
			self.send(self.kafka_config['simulation_receive'], {
				'event': event.SHUTDOWN_SIMULATION,
				'inner_id': key})

			removed[key] = self.agents.pop(key)

		for key in removing:
			removed[key]['process'].join()

		print('removal complete {}'.format(removing))

		return removed

	def receive(self, topic, message):
		"""
		This agent receives two events: ADD_AGENT and REMOVE_AGENT, who's message contents
		match the arguments to the functions above.
		"""

		print('--> {}: {}'.format(topic, message))

		if message['event'] == event.ADD_AGENT:
			self.add_agent(
				message['agent_id'],
				message['agent_type'],
				message['agent_config'])

		elif message['event'] == event.REMOVE_AGENT:
			self.remove_agent(message['agent_prefix'])
