#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

from wholecell.utils.fitting import calcProteinCounts

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

		aaIDs = kb.moleculeGroups.aaIDs
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedAA_IDs]

		## Find the average protein composition and initial number of polymerized AAs

		self.cellCycleLen = kb.doubling_time.asNumber(units.s)

		proteinComposition = kb.process.translation.monomerData["aaCounts"].asNumber()

		initialDryMass = kb.mass.avgCellDryMassInit

		proteinMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].asNumber(units.min) == 60.0
			]["proteinMassFraction"]

		initialProteinMass = initialDryMass * proteinMassFraction

		initialProteinCounts = calcProteinCounts(kb, initialProteinMass)

		initialProteinTranslationRate = initialProteinCounts * (
			np.log(2) / kb.doubling_time + kb.process.translation.monomerData["degRate"]
			).asNumber(1 / units.s)

		initialPolymerizing = np.dot(proteinComposition.T, initialProteinTranslationRate)

		self.initialPolymerizingTotal = initialPolymerizing.sum()

		self.monomerComposition = initialPolymerizing / initialPolymerizing.sum()

		## Energy costs

		self.gtpPerElongation = kb.constants.gtpPerTranslation

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
			self.initialPolymerizingTotal
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
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
