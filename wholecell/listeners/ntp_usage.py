#!/usr/bin/env python

"""
NtpUsage

NtpUsage listener. Tracks NTP usages.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/21/2014
"""

from __future__ import division

import numpy as np
import tables

import wholecell.listeners.listener
from wholecell.reconstruction.fitter import normalize

class NtpUsage(wholecell.listeners.listener.Listener):
	""" NtpUsage """

	_name = 'NtpUsage'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(NtpUsage, self).__init__(*args, **kwargs)

		self.usageUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(NtpUsage, self).initialize(sim, kb)

		self.bulkMolecules = sim.states['BulkMolecules']

		self.states = sim.states

		self.sim = sim

		self.metaboliteIds = ["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"]
		self.metaboliteIdxs = [
			np.where(
				x == self.bulkMolecules._moleculeIDs
				)[0][0] for x in self.metaboliteIds
		]

		self.transcriptionProcessIdx = sim.processes.keys().index(
			"UniqueTranscriptElongation"
			)

		biomassMetIdxs = [
		np.where(
			kb.wildtypeBiomass["metaboliteId"] == x
			)[0][0] for x in self.metaboliteIds
		]

		self.relativeNtpProductionBiomass = normalize(
			kb.wildtypeBiomass["biomassFlux"][biomassMetIdxs].magnitude
			)

	# Allocate memory
	def allocate(self):
		super(NtpUsage, self).allocate()

		self.transcriptionNtpUsageCurrent = np.zeros(
			len(self.metaboliteIds), dtype = np.int64
			)
		self.transcriptionNtpUsageCumulative = np.zeros_like(
			self.transcriptionNtpUsageCurrent
			)

	def update(self):
		self.transcriptionNtpUsageCurrent = (
			self.bulkMolecules._countsAllocatedInitial[
				self.metaboliteIdxs, self.transcriptionProcessIdx
				]
			- self.bulkMolecules._countsAllocatedFinal[
				self.metaboliteIdxs, self.transcriptionProcessIdx
				]
			)

		self.transcriptionNtpUsageCumulative += (
			self.transcriptionNtpUsageCurrent
			)

	def pytablesCreate(self, h5file, expectedRows):

		shape = self.transcriptionNtpUsageCurrent.shape
		# Columns
		d = {
			"time": tables.Int64Col(),
			"transcriptionNtpUsageCurrent": tables.UInt64Col(shape),
			"transcriptionNtpUsageCumulative": tables.UInt64Col(shape),
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
		t.attrs.transcriptionNtpUsageCurrent_units = self.usageUnits
		t.attrs.transcriptionNtpUsageCumulative_units = self.usageUnits
		t.attrs.metaboliteIds = self.metaboliteIds
		t.attrs.relativeNtpProductionBiomass = self.relativeNtpProductionBiomass


	def pytablesAppend(self, h5file):

		t = h5file.get_node("/", self._name)
		entry = t.row

		entry["time"] = self.timeStep()
		entry["transcriptionNtpUsageCurrent"] = self.transcriptionNtpUsageCurrent
		entry["transcriptionNtpUsageCumulative"] = self.transcriptionNtpUsageCumulative

		entry.append()

		t.flush()