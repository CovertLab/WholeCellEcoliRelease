from __future__ import absolute_import, division, print_function

from agent.inner import CellSimulation

DEFAULT_COLOR = [color/255 for color in [0, 128, 255]]


class TransportComposite(CellSimulation):
	'''
	This CellSimulation class is a composite agent, which coordinates the function of multiple sub-agents.

	Its primary function is to initialize the individual sub-agent, adding additional arguments to synchronize_config,
	and then to manage the message passing between them.

	Args:
		network_config (dict) includes:
			* initialize -- a dict mapping each subprocess id to a function for initializing that subprocess
			* message_connections -- a dict for crossing messages with:
					{source_process.source_message: target_process.target_message}
			* function_connections -- a dict for assigning a subprocess function to the composite with:
					{function_name: source_process}
	'''

	def __init__(self, boot_config, synchronize_config, network_config):
		initialize_processes = network_config['initialize']
		function_connections = network_config['function_connections']
		self.connections = network_config['message_connections']

		# initialize message passing between subprocesses
		self.cross_update = {subprocess: {} for subprocess in initialize_processes.keys()}

		# initialize transport
		self.processes = {}
		self.processes['transport'] = initialize_processes['transport'](boot_config, synchronize_config)

		# update ecoli's synchronize_config using transport_reactions_ids from transport
		ecoli_synchronize_config = boot_config
		ecoli_synchronize_config['boundary_reactions'] = self.processes['transport'].transport_reactions_ids

		# initialize ecoli
		self.processes['ecoli'] = initialize_processes['ecoli'](boot_config, ecoli_synchronize_config)

		# use source process' functions as composite function
		# TODO (eran) what happens when one process declares division? Can it divide the rest? How do they share inheritance, etc.
		self.time = self.processes[function_connections['time']].time
		self.divide = self.processes[function_connections['divide']].divide


	def generate_inner_update(self):
		'''
		Sends environment a dictionary with relevant state changes.
		Creates dictionary with messages between processes.
		'''

		# get updates from all of the individual processes
		process_updates = {}
		for process_id, process in self.processes.iteritems():
			process_update = process.generate_inner_update()
			process_updates[process_id] = process_update

		# add to inner update by pulling from individual process updates according to connections
		inner_update = {}
		for source, target in self.connections.iteritems():
			source_process, source_message = source.split('.')
			target_process, target_message = target.split('.')

			if target_process == 'composite':
				inner_update[target_message] = process_updates[source_process][source_message]
			else:
				self.cross_update[target_process][target_message] = process_updates[source_process][source_message]

		# inner updates from composite
		inner_update['color'] = DEFAULT_COLOR

		return inner_update


	def apply_outer_update(self, outer_update):
		for process_id, process in self.processes.iteritems():
			cross_update = self.cross_update[process_id]
			# merge dicts
			process_update = dict(outer_update, **cross_update)  # Python 3 uses dict(**a, **b)
			process.apply_outer_update(process_update)

	def run_incremental(self, run_until):
		for process in self.processes.itervalues():
			process.run_incremental(run_until)

	def finalize(self):
		for process in self.processes.itervalues():
			process.finalize()
