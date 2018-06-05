#!/usr/bin/env python

"""
ReplicationData

Replication fork position listener. Represents position of replication forks over time.

@author: Nick Ruggero
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

		self.sequenceIdx = np.zeros(75, np.int64)
		self.sequenceIdx.fill(PLACE_HOLDER)
		self.sequenceLength = np.zeros(75, np.float64)
		self.sequenceLength.fill(PLACE_HOLDER)

		self.numberOfOric = np.nan
		self.criticalMassPerOriC = 0.
		self.criticalInitiationMass = 0.


	def update(self):
		dnaPolymerases = self.uniqueMolecules.container.objectsInCollection('dnaPolymerase')
		oriCs = self.uniqueMolecules.container.objectsInCollection('originOfReplication')

		self.numberOfOric = len(oriCs)

		if len(dnaPolymerases) > 0:
			sequenceIdx, sequenceLength = dnaPolymerases.attrs(
				"sequenceIdx", "sequenceLength"
				)
			self.sequenceIdx[:] = PLACE_HOLDER
			self.sequenceIdx[:sequenceIdx.size] = sequenceIdx
			self.sequenceLength[:] = np.NAN
			self.sequenceLength[:sequenceLength.size] = sequenceLength
		elif len(dnaPolymerases) == 0:
			self.sequenceIdx[:] = PLACE_HOLDER
			self.sequenceLength[:] = PLACE_HOLDER

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			sequenceIdx = self.sequenceIdx,
			sequenceLength = self.sequenceLength,
			numberOfOric = self.numberOfOric,
			criticalMassPerOriC = self.criticalMassPerOriC,
			criticalInitiationMass = self.criticalInitiationMass,
			)
