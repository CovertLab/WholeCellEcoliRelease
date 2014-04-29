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
		ntpCounts = self.ntps.counts()

		activeRnaPolys = self.activeRnaPolys.molecules()

		if len(activeRnaPolys) == 0:
			return

		assignedNts, requiredNts, massDiffRna = activeRnaPolys.attrs(
			'assignedAUCG', 'requiredAUCG', 'massDiffRna'
			)

		deficitNts = requiredNts - assignedNts

		newlyAssignedNts = polymerizePooledMonomers(
			ntpCounts,
			deficitNts,
			self.elngRate,
			self.randStream,
			useIntegerLinearProgramming = False
			)

		ntpsUsed = newlyAssignedNts.sum(0)

		updatedNts = newlyAssignedNts + assignedNts

		didInitialize = (assignedNts.sum(1) == 0) & (updatedNts.sum(1) > 0)

		updatedMass = massDiffRna + np.dot(newlyAssignedNts, self.ntWeights)

		updatedMass[didInitialize] += self.hydroxylWeight

		activeRnaPolys.attrIs(
			assignedAUCG = updatedNts,
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

		# HAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAACK

		totalRnapSubunits = self.rnapSubunits.total()/np.array([2, 1, 1, 1], np.int) + len(activeRnaPolys)

		totalRnap = np.min(totalRnapSubunits)

		self.rnapSubunits.countsInc(np.fmax(
			((440*np.exp(np.log(2)/3600*self.time()) - totalRnap) * np.array([2, 1, 1, 1], np.int)).astype(np.int),
			0
			))


# TODO: move out of this file
from cvxopt.glpk import ilp, lp
from cvxopt import matrix, sparse, spmatrix

# TODO: line-profile
def polymerizePooledMonomers(monomerCounts, monomerDeficits, maxElongation,
		randStream, useIntegerLinearProgramming = True):

	deficits = monomerDeficits.reshape(-1)

	nMonomers = monomerDeficits.shape[1]
	nPolymers = monomerDeficits.shape[0]

	nNodes = nMonomers + nPolymers # i.e. rows
	nEdges = nMonomers * (1 + nPolymers) + nPolymers # i.e. columns

	nNonzeroEntries = nMonomers + 2*nMonomers*nPolymers + nPolymers

	values = np.empty(nNonzeroEntries)
	rowIndexes = np.empty(nNonzeroEntries, np.int)
	colIndexes = np.empty(nNonzeroEntries, np.int)

	values[:nMonomers] = 1
	values[nMonomers:nMonomers*(1+nPolymers)] = -1
	values[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = 1
	values[-nPolymers:] = -1

	rowIndexes[:nMonomers*(1+nPolymers)] = np.tile(np.arange(nMonomers), 1+nPolymers)
	rowIndexes[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = nMonomers + np.repeat(np.arange(nPolymers), nMonomers)
	rowIndexes[-nPolymers:] = nMonomers + np.arange(nPolymers)

	colIndexes[:nMonomers*(1+nPolymers)] = np.arange(nMonomers*(1+nPolymers))
	colIndexes[nMonomers*(1+nPolymers):nMonomers*(1+2*nPolymers)] = nMonomers + np.arange(nMonomers*nPolymers)
	colIndexes[-nPolymers:] = nMonomers*(1+nPolymers) + np.arange(nPolymers)

	lowerBounds = np.zeros(nEdges)

	upperBounds = np.empty(nEdges)

	upperBounds[:nMonomers] = monomerCounts
	upperBounds[nMonomers:nMonomers+nMonomers*nPolymers] = deficits
	upperBounds[nMonomers+nMonomers*nPolymers:] = maxElongation

	# Add some noise to prevent systematic bias in polymer NT assignment
	# NOTE: still not a great solution, see if a flexFBA solution exists
	objectiveNoise = 0.01 * (randStream.rand(nPolymers) / nPolymers)

	objective = np.zeros(nEdges)
	objective[nMonomers+nMonomers*nPolymers:] = -1 * (1 + objectiveNoise)

	# Setting up (I)LP algorithm:
	# Choose v such that
	# b = Av
	# v < Gh
	# min f^T v
	# (all v integer-valued)

	A = spmatrix(values, rowIndexes, colIndexes)
	b = matrix(np.zeros(nNodes))

	# TODO: build this matrix smarter
	G = spmatrix(
		[1]*nEdges + [-1]*nEdges,
		np.arange(2*nEdges),
		np.tile(np.arange(nEdges), 2)
		)

	h = matrix(np.concatenate([upperBounds, -lowerBounds], axis = 0))

	f = matrix(objective)

	if useIntegerLinearProgramming:
		I = set(xrange(nEdges))
		solution = ilp(f, G, h, A, b, I)

	else:
		solution = lp(f, G, h, A, b)

	fluxes = np.array(solution[1]).astype(np.int).flatten()

	assignments = fluxes[nMonomers:nMonomers+nMonomers*nPolymers]

	monomerAssignments = assignments.reshape(nPolymers, nMonomers)

	return monomerAssignments
