from __future__ import absolute_import, division, print_function

import os
import time
import uuid
import errno
import argparse

from agent.control import AgentControl, AgentCommand

class ShepherdControl(AgentControl):

	"""
	Send messages to the other agents in the system to trigger execution and/or shutdown
	the Outer agent (which sends messages to shutdown all the associated Inner agents) or
	shutdown specific Inner agents directly (which then report back to the Outer agent and
	then terminate).
	"""

	def __init__(self, agent_config):
		super(ShepherdControl, self).__init__(str(uuid.uuid1()), agent_config)

	def add_cell(self, agent_type, agent_config):
		self.add_agent(
			str(uuid.uuid1()),
			agent_type,
			agent_config)

	def lattice_experiment(self, args):
		lattice_id = str(uuid.uuid1())
		self.add_agent(lattice_id, 'lattice', {})

		for index in range(args['number']):
			self.add_cell(args['type'] or 'ecoli', {
				'outer_id': lattice_id,
				'working_dir': args['working_dir']})

	def chemotaxis_experiment(self, args):
		lattice_id = str(uuid.uuid1())
		chemotaxis_config = {
			'run_for' : 1.0,
			'static_concentrations': True,
			'gradient': {'seed': True},
			'diffusion': 0.0,
			'translation_jitter': 0.0,
			'rotation_jitter': 0.0,
			'edge_length': 20.0,
			'patches_per_edge': 30}
		self.add_agent(lattice_id, 'lattice', chemotaxis_config)

		for index in range(args['number']):
			self.add_cell(args['type'] or 'chemotaxis', {
				'outer_id': lattice_id})


class EnvironmentCommand(AgentCommand):
	"""
	Extend `AgentCommand` with new commands related to the lattice and ecoli experiments
	"""

	def __init__(self):
		choices = ['ecoli', 'chemotaxis', 'lattice', 'chemotaxis-experiment']
		description = '''
		Run an agent for the environmental context simulation.
		The commands are:
		`add --id OUTER_ID [--type T] [--variant V] [--index I] [--seed S]` ask the Shepherd to add an agent of type T,
		`experiment [--number N] [--type T] [--media M]` ask the Shepherd to run a lattice environment with N agents of type T in media condition M,
		`pause --id OUTER_ID` pause the simulation,
		`remove --id OUTER_ID` ask all Shepherds to remove agents "ID*",
		`shepherd [--working-dir D]` run a Shepherd agent in this process,
		`shutdown --id OUTER_ID` shut down the environment agent and all its connected agents,
		`trigger --id OUTER_ID` start running the simulation'''

		super(EnvironmentCommand, self).__init__(choices, description)

	def experiment(self, args):
		control = ShepherdControl({'kafka_config': self.kafka_config})
		control.lattice_experiment(args)
		control.shutdown()

	def chemotaxis_experiment(self, args):
		control = ShepherdControl({'kafka_config': self.kafka_config})
		control.chemotaxis_experiment(args)
		control.shutdown()


if __name__ == '__main__':
	command = EnvironmentCommand()
	command.execute()
