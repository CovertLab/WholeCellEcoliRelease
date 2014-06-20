#!/usr/bin/env python

"""
BiomassInternment

A process that only performs requests, interning the molecules over time steps.
Useful for processes that are sensitive to the metabolic state of the cell.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/12/2014
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_INTERN

class BiomassInternment(wholecell.processes.process.Process):
	""" BiomassInternment """

	_name = "BiomassInternment"

	# Constructor
	def __init__(self):
		super(BiomassInternment, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(BiomassInternment, self).initialize(sim, kb)

		moleculeIds = (
			list(kb.cellGlycogenFractionData["metaboliteId"]) +
			list(kb.cellMureinFractionData["metaboliteId"]) +
			list(kb.cellLPSFractionData["metaboliteId"]) +
			list(kb.cellLipidFractionData["metaboliteId"]) +
			list(kb.cellInorganicIonFractionData["metaboliteId"]) +
			list(kb.cellSolublePoolFractionData["metaboliteId"]) + 
			["H2O[c]", "H[c]"]
			)

		self.molecules = self.bulkMoleculesView(moleculeIds)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_INTERN)


	def calculateRequest(self):
		self.molecules.requestAll()


	def evolveState(self):
		pass
