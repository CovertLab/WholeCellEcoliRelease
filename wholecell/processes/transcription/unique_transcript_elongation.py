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
from wholecell.utils.polymerize_pooled_monomers import polymerizePooledMonomers

# TODO: vectorize operations
# TODO: write ILP solver
# TODO: refactor mass calculations
# TODO: confirm reaction stoich
# TODO: resolve mounting process namespace issues
# TODO: use more intelligent requests (enzyme/metabolite-limited)
# TODO: resolve nucleotide pooling (fix biomass function?)

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

		# TODO: refactor mass updates

		# TOKB
		self.ntWeights = np.array([
			345.20, # A
			322.17, # U
			321.18, # C
			361.20  # G
			]) - 17.01 # weight of a hydroxyl

		# TOKB
		self.hydroxylWeight = 17.01 # counted once for the end of the polymer

		self.ntWeights *= 1e15/6.022e23
		self.hydroxylWeight *= 1e15/6.022e23

		# Views

		self.activeRnaPolys = self.uniqueMoleculesView('activeRnaPoly')
		self.bulkRnas = self.bulkMoleculesView(self.rnaIds)

		self.ntps = self.bulkMoleculesView(["ATP[c]", "CTP[c]", "GTP[c]", "UTP[c]"])
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
		ntpCounts = self.ntps.counts()

		activeRnaPolys = self.activeRnaPolys.molecules()

		if len(activeRnaPolys) == 0:
			return

		assignedNts, requiredNts, massDiffRna = activeRnaPolys.attrs(
			'assignedACGU', 'requiredACGU', 'massDiffRna'
			)

		deficitNts = requiredNts - assignedNts

		updatedNts = assignedNts.copy()
		ntpsUsed = np.zeros_like(ntpCounts)

		wholecell.utils.polymerize.polymerize(
			self.elngRate, deficitNts, requiredNts, ntpCounts,
			updatedNts, ntpsUsed, self.seed
			)

		assert np.all(updatedNts <= requiredNts), "Transcripts got elongated more than possible!"

		# newlyAssignedNts = polymerizePooledMonomers(
		# 	ntpCounts,
		# 	deficitNts,
		# 	self.elngRate,
		# 	self.randStream,
		# 	useIntegerLinearProgramming = False
		# 	)

		# ntpsUsed = newlyAssignedNts.sum(axis = 0)

		# updatedNts = newlyAssignedNts + assignedNts

		didInitialize = (
			(assignedNts.sum(axis = 1) == 0) &
			(updatedNts.sum(axis = 1) > 0)
			)

		updatedMass = massDiffRna + np.dot(
			(updatedNts - assignedNts), self.ntWeights
			)

		updatedMass[didInitialize] += self.hydroxylWeight

		activeRnaPolys.attrIs(
			assignedACGU = updatedNts,
			massDiffRna = updatedMass
			)

		terminatedRnas = np.zeros_like(self.bulkRnas.counts())

		didTerminate = (requiredNts == updatedNts).all(axis = 1)

		for moleculeIndex, molecule in enumerate(activeRnaPolys):
			if didTerminate[moleculeIndex]:
				terminatedRnas[molecule.attr('rnaIndex')] += 1
				self.activeRnaPolys.moleculeDel(molecule)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = ntpsUsed.sum()

		self.ntps.countsDec(ntpsUsed)

		self.bulkRnas.countsIs(terminatedRnas)

		self.rnapSubunits.countsInc(
			nTerminated * np.array([2, 1, 1, 1], np.int) # complex subunit stoich
			)

		self.h2o.countDec(nInitialized)
		self.proton.countInc(nInitialized)

		self.ppi.countInc(nElongations)
