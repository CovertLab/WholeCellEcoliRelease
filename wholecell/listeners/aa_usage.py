#!/usr/bin/env python

"""
AAUsage

AAUsage listener. Tracks amino acid usages.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener
from wholecell.reconstruction.fitter import normalize

class AAUsage(wholecell.listeners.listener.Listener):
	""" AAUsage """

	_name = 'AAUsage'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(AAUsage, self).__init__(*args, **kwargs)

		self.usageUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(AAUsage, self).initialize(sim, kb)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.states = sim.states

		self.sim = sim

		self.metaboliteIds = kb.aaIDs[:]
		self.metaboliteIdxs = [ # TODO: use a bulk container view?
			np.where(
				x == self.bulkMolecules._moleculeIDs
				)[0][0] for x in self.metaboliteIds
			]

		self.translationProcessIdx = sim.processes.keys().index(
			"UniquePolypeptideElongation"
			)

		biomassMetIdxs = [
		np.where(
			kb.wildtypeBiomass["metaboliteId"] == x
			)[0][0] for x in self.metaboliteIds
		]

		self.relativeAAProductionBiomass = normalize(
			kb.wildtypeBiomass["biomassFlux"][biomassMetIdxs].magnitude
			)

		self.relativeAaUsage = normalize(np.dot(
			kb.monomerData["aaCounts"].T,
			kb.rnaExpression["expression"][kb.rnaIndexToMonomerMapping]
			))

	# Allocate memory
	def allocate(self):
		super(AAUsage, self).allocate()

		self.translationAAUsageCurrent = np.zeros(
			len(self.metaboliteIds), dtype = np.int64
			)
		self.translationAAUsageCumulative = np.zeros_like(
			self.translationAAUsageCurrent
			)

	def update(self):
		self.translationAAUsageCurrent = (
			self.bulkMolecules._countsAllocatedInitial[
				self.metaboliteIdxs, self.translationProcessIdx
				]
			- self.bulkMolecules._countsAllocatedFinal[
				self.metaboliteIdxs, self.translationProcessIdx
				]
			)

		self.translationAAUsageCumulative += (
			self.translationAAUsageCurrent
			)

	def pytablesCreate(self, h5file, expectedRows):

		shape = self.translationAAUsageCurrent.shape
		# Columns
		d = {
			"time": tables.Int64Col(),
			"translationAAUsageCurrent": tables.UInt64Col(shape),
			"translationAAUsageCumulative": tables.UInt64Col(shape),
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
		t.attrs.translationAAUsageCurrent_units = self.usageUnits
		t.attrs.translationAAUsageCumulative_units = self.usageUnits
		t.attrs.metaboliteIds = self.metaboliteIds
		t.attrs.relativeAAProductionBiomass = self.relativeAAProductionBiomass
		t.attrs.relativeAaUsage = self.relativeAaUsage


	def pytablesAppend(self, h5file):

		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.timeStep()
		entry["translationAAUsageCurrent"] = self.translationAAUsageCurrent
		entry["translationAAUsageCumulative"] = self.translationAAUsageCumulative

		entry.append()

		t.flush()