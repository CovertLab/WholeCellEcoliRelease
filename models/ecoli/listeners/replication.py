#!/usr/bin/env python

"""
Replication

Replication fork position listener. Represents position of replication forks over time.

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/13/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class Replication(wholecell.listeners.listener.Listener):
	""" Replication """

	_name = 'Replication'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(Replication, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(Replication, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states['UniqueMolecules']

		self.dnaPolyData = None


	# Allocate memory
	def allocate(self):
		super(Replication, self).allocate()


	def update(self):
		dnaPolymerases = self.uniqueMolecules.container.objectsInCollection('dnaPolymerase')

		if len(dnaPolymerases) > 0:
			self.dnaPolyData = dnaPolymerases.attrsAsStructArray(
				"_uniqueId",
				"sequenceIdx",
				"sequenceLength",
				"replicationRound",
				"replicationDivision"
				)
		else:
			# TODO: John rewrite this attrsAsStructArray function to build the struct array with
			# cached data and then populate it in a second functin. Then we can just save the 
			# empty one here.
			self.dnaPolyData = np.zeros(0, dtype = [('_uniqueId', 'S40'), ('sequenceIdx', '<i8'), ('sequenceLength', '?')])

	def tableCreate(self, tableWriter):
		pass


	def tableAppend(self, tableWriter):
		tableWriter.append(
			dnaPolyData = self.dnaPolyData
			)
