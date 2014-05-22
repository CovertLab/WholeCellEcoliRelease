#!/usr/bin/env python

"""
RibosomeStalling

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/14
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener

from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class RibosomeStalling(wholecell.listeners.listener.Listener):
	""" RibosomeStalling """

	_name = 'RibosomeStalling'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(RibosomeStalling, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(RibosomeStalling, self).initialize(sim, kb)

		self.uniqueMolecules = sim.states["UniqueMolecules"]

		self.elngRate = 16 # TODO: get from KB


	# Allocate memory
	def allocate(self):
		super(RibosomeStalling, self).allocate()

		self.stalledRibosomes = np.zeros(
				0,
				dtype = [("t", "i8"), ("stalls", "i8"), ("proteinIndex", "i8")]
				)


	def initialUpdate(self):
		self._recordAssignments()


	def update(self):
		if len(self.ribosomes):
			# Find difference in elongation
			updatedAAs = self.ribosomes.attr("assignedAAs").sum(1)

			actualElongations = updatedAAs - self.assignedAAs

			stalls = self.expectedElongations - actualElongations

			proteinIndexes = self.ribosomes.attr("proteinIndex")

			self.stalledRibosomes = np.zeros(
				stalls.size,
				dtype = self.stalledRibosomes.dtype
				)

			self.stalledRibosomes["t"] = self.timeStep()
			self.stalledRibosomes["stalls"] = stalls
			self.stalledRibosomes["proteinIndex"] = proteinIndexes

			if VERBOSE:
				print (stalls > 0).sum(), stalls.mean(), stalls.std(), stalls.sum()

		self._recordAssignments()


	def _recordAssignments(self):
		# Record new quantities for next simulation step
		self.ribosomes = self.uniqueMolecules.container.objectsInCollection("activeRibosome")

		if len(self.ribosomes) == 0:
			return

		self.assignedAAs = self.ribosomes.attr("assignedAAs").sum(1)
		self.expectedElongations = np.fmin(
			self.ribosomes.attr("requiredAAs").sum(1) - self.assignedAAs,
			self.elngRate
			)


	def pytablesCreate(self, h5file, expectedRows):
		table = h5file.create_table(
			h5file.root,
			self._name,
			self.stalledRibosomes.dtype,
			title = self._name,
			filters = tables.Filters(complevel = 9, complib="zlib"),
			)


	def pytablesAppend(self, h5file):
		table = h5file.get_node("/", self._name)

		table.append(self.stalledRibosomes)

		table.flush()
