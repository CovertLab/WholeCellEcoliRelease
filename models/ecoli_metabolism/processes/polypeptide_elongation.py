#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units
from reconstruction.ecoli.fitter import calcProteinCounts

class PolypeptideElongation(wholecell.processes.process.Process):
	""" PolypeptideElongation """

	_name = "PolypeptideElongation"

	# Constructor
	def __init__(self):
		super(PolypeptideElongation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(PolypeptideElongation, self).initialize(sim, kb)

		# TODO: account for the polymer ends

		# Load constants

		## Molecule IDs

		aaIDs = kb.aaIDs
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedAA_IDs]

		## Find the average protein composition and initial number of polymerized AAs

		self.cellCycleLen = kb.cellCycleLen.asUnit(units.s).asNumber()

		proteinComposition = kb.monomerData["aaCounts"].asNumber()

		initialDryMass = kb.avgCellDryMassInit

		proteinMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
			]["proteinMassFraction"]

		initialProteinMass = initialDryMass * proteinMassFraction

		proteinCounts = calcProteinCounts(kb, initialProteinMass)

		polymerizedCounts = np.dot(proteinComposition.T, proteinCounts)

		self.monomerComposition = polymerizedCounts / polymerizedCounts.sum()

		self.initialAverageMonomerCounts = polymerizedCounts.sum()

		## Energy costs

		self.gtpPerElongation = kb.gtpPerTranslation

		# Create views on state

		self.polymerized = self.bulkMoleculesView(polymerizedIDs)

		self.aas = self.bulkMoleculesView(aaIDs)
		self.h2o = self.bulkMoleculeView("H2O[c]")
		self.gtp = self.bulkMoleculeView("GTP[c]")

		self.gdp = self.bulkMoleculeView("GDP[c]")
		self.pi = self.bulkMoleculeView("PI[c]")
		self.h = self.bulkMoleculeView("H[c]")


	def calculateRequest(self):
		totalMonomers = np.int64(stochasticRound(
			self.randomState,
			self.initialAverageMonomerCounts
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
			* (np.exp(np.log(2) / self.cellCycleLen * self.timeStepSec) - 1)
			))

		aasRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		gtpRequested = np.int64(stochasticRound(
			self.randomState,
			self.gtpPerElongation * aasRequested.sum()
			))

		h2oRequested = gtpRequested - aasRequested.sum()

		if (
			(aasRequested > self.aas.total()).any()
			or gtpRequested > self.gtp.total()
			or h2oRequested > self.h2o.total()
				):

			# TODO: flag simulation instead of printing
			print "{} is metabolically limited".format(self.name())

		self.aas.requestIs(aasRequested)
		self.gtp.requestIs(gtpRequested)
		self.h2o.requestIs(h2oRequested)


	def evolveState(self):
		aaCounts = self.aas.counts()

		self.polymerized.countsInc(aaCounts)

		self.h2o.countInc(aaCounts.sum())

		gtpCount = self.gtp.count()

		self.gdp.countInc(gtpCount)
		self.pi.countInc(gtpCount)
		self.h.countInc(gtpCount)

		self.h2o.countDec(gtpCount)

		self.aas.countsIs(0)
		self.gtp.countIs(0)
