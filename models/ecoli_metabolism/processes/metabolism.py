#!/usr/bin/env python

from __future__ import division

from itertools import izip

import numpy as np

import wholecell.processes.process

from wholecell.utils.random import stochasticRound
from wholecell.utils.constants import REQUEST_PRIORITY_METABOLISM
from wholecell.utils import units

class Metabolism(wholecell.processes.process.Process):
	""" Metabolism """

	_name = "Metabolism"

	# Constructor
	def __init__(self):
		super(Metabolism, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(Metabolism, self).initialize(sim, kb)

		# Load constants

		self.nAvogadro = kb.constants.nAvogadro.asNumber(1 / units.mol)
		self.cellDensity = kb.constants.cellDensity.asNumber(units.g/units.L)
		
		self.metabolitePoolIDs = kb.process.metabolism.metabolitePoolIDs
		self.targetConcentrations = kb.process.metabolism.metabolitePoolConcentrations.asNumber(units.mol/units.L)
		
		# Create views on state

		self.metabolites = self.bulkMoleculesView(self.metabolitePoolIDs)

		# Set the priority to a low value

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_METABOLISM)


	def calculateRequest(self):
		self.metabolites.requestAll()


	def evolveState(self):
		# Set metabolite pools to target concentrations
		
		metCountsInit = self.metabolites.counts()

		cellMass = self.readFromListener("Mass", "cellMass") * 1e-15

		cellVolume = cellMass / self.cellDensity

		countsToMolar = 1 / (self.nAvogadro * cellVolume)

		metConcInit = metCountsInit * countsToMolar

		metCountsFinal = self.targetConcentrations / countsToMolar

		self.metabolites.countsIs(np.int64(stochasticRound(
			self.randomState,
			metCountsFinal
			)))

		self.writeToListener("ConcentrationChange", "concentrationChange", 
			self.targetConcentrations - metConcInit)
