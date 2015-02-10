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

import wholecell.listeners.listener
from wholecell.utils.fitting import normalize
from wholecell.utils import units

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
			"PolypeptideElongation"
			)

		biomassMetIdxs = [
		np.where(
			kb.wildtypeBiomass["metaboliteId"] == x
			)[0][0] for x in self.metaboliteIds
		]

		self.relativeAAProductionBiomass = normalize(
			kb.wildtypeBiomass["biomassFlux"][biomassMetIdxs].asNumber()
			)

		self.relativeAaUsage = normalize(np.dot(
			kb.monomerData["aaCounts"].asNumber().T,
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

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			translationAAUsageCurrent_units = self.usageUnits,
			translationAAUsageCumulative_units = self.usageUnits,
			metaboliteIds = self.metaboliteIds,
			relativeAAProductionBiomass = self.relativeAAProductionBiomass.tolist(),
			relativeAaUsage = self.relativeAaUsage.tolist()
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			translationAAUsageCurrent = self.translationAAUsageCurrent,
			translationAAUsageCumulative = self.translationAAUsageCumulative,
			)
