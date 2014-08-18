#!/usr/bin/env python

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.random import stochasticRound
from wholecell.utils import units

from reconstruction.ecoli.fitter import countsFromMassAndExpression

class TranscriptElongation(wholecell.processes.process.Process):
	""" TranscriptElongation """

	_name = "TranscriptElongation"

	def __init__(self):
		super(TranscriptElongation, self).__init__()


	def initialize(self, sim, kb):
		super(TranscriptElongation, self).initialize(sim, kb)

		# TODO: account for the polymer ends

		# Load constants

		## Molecule IDs

		ntpIDs = kb.ntpIds
		polymerizedIDs = [id_ + "[c]" for id_ in kb.polymerizedNT_IDs]

		## Find the magnitude and composition of transcription

		rnaComposition = kb.rnaData["countsACGU"].asNumber()

		initialDryMass = kb.avgCellDryMassInit

		rnaMassFraction = kb.cellDryMassComposition[
			kb.cellDryMassComposition["doublingTime"].asUnit(units.min).asNumber() == 60.0
			]["rnaMassFraction"]

		initialRnaMass = initialDryMass * rnaMassFraction

		initialRnaCounts = countsFromMassAndExpression(
			initialRnaMass.asNumber(units.g),
			kb.rnaData["mw"].asNumber(units.g / units.mol),
			kb.rnaExpression["expression"],
			kb.nAvogadro.asNumber(1 / units.mol)
			) * kb.rnaExpression["expression"]

		initialRnaTranscriptionRate = initialRnaCounts * (
			np.log(2) / kb.cellCycleLen + kb.rnaData["degRate"]
			).asNumber(1 / units.s)

		initialPolymerizing = np.dot(rnaComposition.T, initialRnaTranscriptionRate)

		self.cellCycleLen = kb.cellCycleLen.asNumber(units.s)

		self.initialPolymerizingTotal = initialPolymerizing.sum()

		self.monomerComposition = initialPolymerizing / initialPolymerizing.sum()

		# Create views on state

		self.polymerized = self.bulkMoleculesView(polymerizedIDs)

		self.ntps = self.bulkMoleculesView(ntpIDs)
		self.ppi = self.bulkMoleculeView("PPI[c]")


	def calculateRequest(self):
		totalMonomers = np.int64(stochasticRound(
			self.randomState,
			self.initialPolymerizingTotal
			* np.exp(np.log(2) / self.cellCycleLen * self.time())
			))

		ntpsRequested = self.randomState.multinomial(
			totalMonomers,
			self.monomerComposition
			)

		if (ntpsRequested > self.ntps.total()).any():
			# TODO: flag simulation instead of printing
			print "{} is metabolically limited".format(self.name())

		self.ntps.requestIs(ntpsRequested)


	def evolveState(self):
		ntpCounts = self.ntps.counts()

		self.polymerized.countsInc(ntpCounts)

		self.ppi.countInc(ntpCounts.sum())

		self.ntps.countsIs(0)
