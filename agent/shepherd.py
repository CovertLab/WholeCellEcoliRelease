from __future__ import absolute_import, division, print_function

import multiprocessing as mp
import agent.event as event
from agent.agent import Agent

class AgentShepherd(Agent):
	def __init__(self, agent_id, kafka_config, agent_initializers):
		self.agents = {}
		self.agent_initializers = agent_initializers

		kafka_config['subscribe_topics'] = [kafka_config['shepherd_control']]
		super(AgentShepherd, self).__init__(agent_id, kafka_config)

	def initialize(self):
		print('agent shepherd waiting')

	def add_agent(self, agent_id, agent_type, agent_config):
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
		removing = filter(lambda key: key.startswith(agent_prefix), self.agents.iterkeys())
		removed = {}
		for key in removing:
			self.send(self.kafka_config['simulation_receive'], {
				'event': event.SHUTDOWN_SIMULATION,
				'inner_id': key})

			removed[key] = self.agents.pop(key)

		return removed

	def receive(self, topic, message):
		print('--> {}: {}'.format(topic, message))

		if message['event'] == event.ADD_AGENT:
			self.add_agent(
				message['agent_id'],
				message['agent_type'],
				message['agent_config'])

		elif message['event'] == event.REMOVE_AGENT:
			self.remove_agent(message['agent_prefix'])
