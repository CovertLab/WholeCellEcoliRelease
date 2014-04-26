#!/usr/bin/env python

"""
UniqueTranscriptElongation

Transcription elongation sub-model.

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/26/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process

class UniqueTranscriptElongation(wholecell.processes.process.Process):
	""" UniqueTranscriptElongation """

	_name = "UniqueTranscriptElongation"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None
		self.rnaIds = None
		self.rnaSynthProb = None

		super(UniqueTranscriptElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniqueTranscriptElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = kb.rnaPolymeraseElongationRate.to('nucleotide / s').magnitude

		enzIds = ["EG10893-MONOMER[c]", "RPOB-MONOMER[c]", "RPOC-MONOMER[c]", "RPOD-MONOMER[c]"]

		self.rnaIds = kb.rnaData['id']

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)

		self.ntps = self.bulkMoleculesView(["ATP[c]", "UTP[c]", "CTP[c]", "GTP[c]"])
		self.ppi = self.bulkMoleculeView('PPI[c]')
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proton = self.bulkMoleculeView('H[c]')

		self.rnapSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.activeRnaPolys.requestAll()

		self.ntps.requestAll()

		self.h2o.requestIs(self.ntps.total().sum()) # this drastically overestimates water assignment


	# Calculate temporal evolution
	def evolveState(self):
		# TODO: implement the following ILP algorithm

		# Assign NTs to each rnaTranscript molecule...
		# ...maximizing the total nucelotides assigned
		# ...constrained to non-negative NT assignments
		# ...constraiend to no more than <elongation rate> total assignments
		# ...constrained to integer values of assignments
		# ...constrained to the maximal NT deficit of each transcript

		activeRnaPolys = list(self.activeRnaPolys.molecules())

		self.randStream.numpyShuffle(activeRnaPolys)

		ntpCounts = self.ntps.counts()
		
		# TODO: vectorize this operation (requires some new unique object accesors)

		terminatedRnas = np.zeros_like(self.bulkRnas.counts())

		freeRnapSubunits = np.zeros_like(self.rnapSubunits.counts())

		nInitialized = 0
		nElongations = 0

		for activeRnaPoly in activeRnaPolys:
			if ntpCounts.sum() == 0:
				break

			assignedAUCG, requiredAUCG = activeRnaPoly.attrs(
				'assignedAUCG',	'requiredAUCG')

			ntDeficit = (requiredAUCG - assignedAUCG).astype(np.float) # for division

			extendedAUCG = np.fmin(
				ntDeficit,
				np.fmin(
					ntDeficit / ntDeficit.sum() * self.elngRate * self.timeStepSec,
					ntpCounts
					)
				).astype(np.int)

			ntpCounts -= extendedAUCG

			updatedAUCG = assignedAUCG + extendedAUCG

			activeRnaPoly.attrIs(
				assignedAUCG = updatedAUCG
				)

			# TODO: update mass

			# TODO: check this elongation reaction stoich
			nElongations += extendedAUCG.sum()

			if assignedAUCG.sum() == 0 and updatedAUCG.sum() > 0:
				nInitialized += 1

			if (updatedAUCG == requiredAUCG).all():
				terminatedRnas[activeRnaPoly.attr('rnaIndex')] += 1

				self.activeRnaPolys.moleculeDel(activeRnaPoly)


		self.ntps.countsIs(ntpCounts)

		self.bulkRnas.countsInc(terminatedRnas)

		self.rnapSubunits.countsInc(
			terminatedRnas.sum() * np.array([2, 1, 1, 1]) # complex subunit stoich
			)

		self.h2o.countDec(nInitialized)
		# self.proton.countInc(nInitialized)

		# self.ppi.countInc(nElongations)

		print nInitialized, nElongations, terminatedRnas.sum()
