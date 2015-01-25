#!/usr/bin/env python

"""
ReplicationForkPosition

Replication fork position listener. Represents position of replication forks over time.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class ReplicationForkPosition(wholecell.listeners.listener.Listener):
	""" ReplicationForkPosition """

	_name = 'ReplicationForkPosition'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ReplicationForkPosition, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(ReplicationForkPosition, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states['UniqueMolecules']

		self.dnaPolyData = None


	# Allocate memory
	def allocate(self):
		super(ReplicationForkPosition, self).allocate()


	def update(self):
		self.dnaPolyData = self.uniqueMolecules.container.objectsInCollection(
			'dnaPolymerase').attrsAsStructArray(
			"_uniqueId",
			"chromosomeLocation",
			"isLeading"
			)


	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			dnaPolyData = self.dnaPolyData
			)
