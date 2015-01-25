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

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

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
			"TranscriptElongation"
			)

		biomassMetIdxs = [
		np.where(
			kb.wildtypeBiomass["metaboliteId"] == x
			)[0][0] for x in self.metaboliteIds
		]

		self.relativeNtpProductionBiomass = normalize(
			kb.wildtypeBiomass["biomassFlux"][biomassMetIdxs].asNumber()
			)

		self.relativeNtpUsage = normalize(
			np.dot(kb.rnaData["countsACGU"].asNumber().T, kb.rnaData["synthProb"])
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

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			transcriptionNtpUsageCurrent_units = self.usageUnits,
			transcriptionNtpUsageCumulative_units = self.usageUnits,
			metaboliteIds = self.metaboliteIds,
			relativeNtpProductionBiomass = self.relativeNtpProductionBiomass,
			relativeNtpUsage = self.relativeNtpUsage,
			)


	def tableAppend(self, tableWriter):
		tableWriter(
			time = self.time(),
			timeStep = self.timeStep(),
			transcriptionNtpUsageCurrent = self.transcriptionNtpUsageCurrent,
			transcriptionNtpUsageCumulative = self.transcriptionNtpUsageCumulative,
			)
