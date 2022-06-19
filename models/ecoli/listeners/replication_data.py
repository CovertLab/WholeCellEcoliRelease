"""
ReplicationData

Replication listener. Records dynamics related to replication.
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class ReplicationData(wholecell.listeners.listener.Listener):
	""" ReplicationData """

	_name = 'ReplicationData'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ReplicationData, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ReplicationData, self).initialize(sim, sim_data)

		self.uniqueMolecules = sim.internal_states['UniqueMolecules']

	# Allocate memory
	def allocate(self):
		super(ReplicationData, self).allocate()

		self.fork_coordinates = np.array([], np.int64)
		self.fork_domains = np.array([], np.int32)
		self.fork_unique_index = np.array([], np.int64)

		self.numberOfOric = np.nan
		self.criticalMassPerOriC = 0.
		self.criticalInitiationMass = 0.

		self.total_DnaA_boxes = 0
		self.free_DnaA_boxes = 0


	def update(self):
		active_replisomes = self.uniqueMolecules.container.objectsInCollection('active_replisome')
		oriCs = self.uniqueMolecules.container.objectsInCollection('oriC')

		self.numberOfOric = len(oriCs)

		self.fork_coordinates, self.fork_domains, self.fork_unique_index = active_replisomes.attrs(
			"coordinates", "domain_index", "unique_index")

		DnaA_boxes = self.uniqueMolecules.container.objectsInCollection('DnaA_box')
		DnaA_box_bound = DnaA_boxes.attrs('DnaA_bound')

		self.total_DnaA_boxes = len(DnaA_boxes)
		self.free_DnaA_boxes = np.count_nonzero(np.logical_not(DnaA_box_bound))


	def tableCreate(self, tableWriter):
		tableWriter.set_variable_length_columns(
			'fork_coordinates',
			'fork_domains',
			'fork_unique_index',
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			fork_coordinates = self.fork_coordinates,
			fork_domains = self.fork_domains,
			fork_unique_index = self.fork_unique_index,
			numberOfOric = self.numberOfOric,
			criticalMassPerOriC = self.criticalMassPerOriC,
			criticalInitiationMass = self.criticalInitiationMass,
			free_DnaA_boxes = self.free_DnaA_boxes,
			total_DnaA_boxes = self.total_DnaA_boxes,
			)
