#!/usr/bin/env python

"""
ReplicationData

Replication fork position listener. Represents position of replication forks
over time.

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

PLACE_HOLDER = -1

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

		self.fork_coordinates = np.full(75, PLACE_HOLDER, np.float64)

		self.numberOfOric = np.nan
		self.criticalMassPerOriC = 0.
		self.criticalInitiationMass = 0.


	def update(self):
		active_replisomes = self.uniqueMolecules.container.objectsInCollection('active_replisome')
		oriCs = self.uniqueMolecules.container.objectsInCollection('originOfReplication')

		self.numberOfOric = len(oriCs)

		if len(active_replisomes) > 0:
			fork_coordinates = active_replisomes.attr("coordinates")
			self.fork_coordinates[:] = PLACE_HOLDER
			self.fork_coordinates[:fork_coordinates.size] = fork_coordinates
		else:
			self.fork_coordinates[:] = PLACE_HOLDER

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			fork_coordinates = self.fork_coordinates,
			numberOfOric = self.numberOfOric,
			criticalMassPerOriC = self.criticalMassPerOriC,
			criticalInitiationMass = self.criticalInitiationMass,
			)
