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
import tables

import wholecell.listeners.listener

class ReplicationForkPosition(wholecell.listeners.listener.Listener):
	""" ReplicationForkPosition """

	_name = 'ReplicationForkPosition'

	# Constructor
	def __init__(self, *args, **kwargs):
		self.positionUnits = 'nt'

		super(ReplicationForkPosition, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(ReplicationForkPosition, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states['UniqueMolecules']

		self.states = sim.states

	# Allocate memory
	def allocate(self):
		super(ReplicationForkPosition, self).allocate()

		self.fork0position = None
		self.fork1position = None

	def _resetPositions(self):
		self.fork0position = 0
		self.fork1position = 0

	def update(self):
		self._resetPositions()
		self.dnaPolymerases = self.uniqueMolecules.container.objectsInCollection('activeDnaPolymerase')

		for i,polymerase in enumerate(self.dnaPolymerases):
			if i == 0:
				self.fork0position = polymerase.attr('chromosomeLocation')
			elif i == 1:
				self.fork1position = polymerase.attr('chromosomeLocation')
			else:
				raise Exception, 'This class is a hack and you have exposed it!\n'


	def pytablesCreate(self, h5file, expectedRows):
		# Columns
		d = {
			"fork0position": tables.Int64Col(),
			"fork1position": tables.Int64Col(),
			}

		# Create table
		# TODO: Add compression options (using filters)
		t = h5file.create_table(
			h5file.root,
			self._name,
			d,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			expectedrows = expectedRows
			)

		# Store units as metadata
		t.attrs.position_units = self.positionUnits

	def pytablesAppend(self, h5file):
		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["fork0position"] = self.fork0position
		entry["fork1position"] = self.fork1position

		entry.append()

		t.flush()
