#!/usr/bin/env python

"""
ReplicationData

Replication listener. Records dynamics related to replication.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

MAX_ACTIVE_REPLISOMES = 32

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

		self.fork_coordinates = np.full(MAX_ACTIVE_REPLISOMES, np.nan, np.float64)
		self.fork_domains = np.full(MAX_ACTIVE_REPLISOMES, np.nan, np.float64)
		self.fork_unique_index = np.full(MAX_ACTIVE_REPLISOMES, np.nan, np.float64)

		self.numberOfOric = np.nan
		self.criticalMassPerOriC = 0.
		self.criticalInitiationMass = 0.

		self.free_DnaA_boxes = 0
		self.total_DnaA_boxes = 0


	def update(self):
		self.fork_coordinates[:] = np.nan
		self.fork_domains[:] = np.nan
		self.fork_unique_index[:] = np.nan

		active_replisomes = self.uniqueMolecules.container.objectsInCollection('active_replisome')
		oriCs = self.uniqueMolecules.container.objectsInCollection('oriC')

		self.numberOfOric = len(oriCs)

		if len(active_replisomes) > 0:
			fork_coordinates, fork_domains, fork_unique_index = active_replisomes.attrs(
				"coordinates", "domain_index", "unique_index"
				)

			self.fork_coordinates[:fork_coordinates.size] = fork_coordinates
			self.fork_domains[:fork_domains.size] = fork_domains
			self.fork_unique_index[:fork_unique_index.size] = fork_unique_index


	def tableCreate(self, tableWriter):
		pass

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
