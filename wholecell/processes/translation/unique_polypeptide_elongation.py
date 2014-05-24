#!/usr/bin/env python

"""
UniquePolypeptideElongation

Translation elongation sub-model.

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 4/30/14
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
import wholecell.utils.polymerize
from wholecell.states.bulk_molecules import calculatePartition

class UniquePolypeptideElongation(wholecell.processes.process.Process):
	""" UniquePolypeptideElongation """

	_name = "UniquePolypeptideElongation"

	# Constructor
	def __init__(self):
		# Constants
		self.elngRate = None
		self.proteinIds = None

		super(UniquePolypeptideElongation, self).__init__()


	# Construct object graph
	def initialize(self, sim, kb):
		super(UniquePolypeptideElongation, self).initialize(sim, kb)

		# Load parameters

		self.elngRate = float(kb.ribosomeElongationRate.to('amino_acid / s').magnitude)

		enzIds = ["RRLA-RRNA[c]", "RRSA-RRNA[c]", "RRFA-RRNA[c]"]

		self.proteinIds = kb.monomerData['id']

		aaIds = kb.aaIDs[:]

		# TODO: Remove hack of deleting selenocysteine this way
		selenocysteineIdx = aaIds.index("SEC-L[c]")
		del aaIds[selenocysteineIdx]

		# # TODO: refactor mass updates

		self.h2oWeight = (
			kb.bulkMolecules[
				kb.bulkMolecules["moleculeId"] == "H2O[c]"
				]["mass"].to("fg / mole").magnitude /
			kb.nAvogadro.to("1 / mole").magnitude
			)

		self.aaWeights = np.array([
			kb.bulkMolecules[
				kb.bulkMolecules["moleculeId"] == x
				]["mass"].to("fg / mole").magnitude /
			kb.nAvogadro.to("1 / mole").magnitude
			for x in kb.aaIDs
			if len(kb.bulkMolecules[kb.bulkMolecules["moleculeId"] == x]["mass"])
			])

		self.aaWeightsIncorporated = self.aaWeights - self.h2oWeight

		# Views

		self.activeRibosomes = self.uniqueMoleculesView('activeRibosome')
		self.bulkMonomers = self.bulkMoleculesView(self.proteinIds)

		self.aas = self.bulkMoleculesView(aaIds)
		self.h2o = self.bulkMoleculeView('H2O[c]')

		self.ribosomeSubunits = self.bulkMoleculesView(enzIds)


	def calculateRequest(self):
		self.activeRibosomes.requestAll()

		activeRibosomes = self.activeRibosomes.allMolecules()

		assignedAAs, requiredAAs = activeRibosomes.attrs("assignedAAs", "requiredAAs")
		deficitAAs = requiredAAs - assignedAAs

		if deficitAAs.size <= 0:
			return

		approxUsage = np.zeros_like(deficitAAs)
		f = deficitAAs / np.tile(deficitAAs.sum(axis = 1).reshape(-1, 1).astype("float64"), (1, 20))

		approxUsage[deficitAAs.sum(axis = 1) <= self.elngRate] = deficitAAs[
			deficitAAs.sum(axis = 1) <= self.elngRate
		]
		approxUsage[deficitAAs.sum(axis = 1) > self.elngRate] = np.ceil(
			f[deficitAAs.sum(axis = 1) > self.elngRate] * self.elngRate
			)

		self.aas.requestIs(
			np.fmin(self.aas.total(), approxUsage.sum(axis = 0))
			)


	# Calculate temporal evolution
	def evolveState(self):
		aaCounts = self.aas.counts()

		activeRibosomes = self.activeRibosomes.molecules()

		if len(activeRibosomes) == 0:
			return

		assignedAAs, requiredAAs, massDiffProtein = activeRibosomes.attrs(
			'assignedAAs', 'requiredAAs', 'massDiffProtein'
			)

		deficitAAs = requiredAAs - assignedAAs

		updatedAAs = assignedAAs.copy()
		aasUsed = np.zeros_like(aaCounts)

		# TODO: perform this in the C function/wrapper so it applies universally
		permutationIndexes = self.randStream.randStream.permutation(assignedAAs.shape[0])
		inversePermutationIndexes = np.argsort(permutationIndexes)

		deficitAAs = deficitAAs[permutationIndexes, :]
		requiredAAs = requiredAAs[permutationIndexes, :]
		updatedAAs = updatedAAs[permutationIndexes, :]

		wholecell.utils.polymerize.polymerize(
			self.elngRate, deficitAAs, requiredAAs, aaCounts,
			updatedAAs, aasUsed, self.seed
			)

		deficitAAs = deficitAAs[inversePermutationIndexes, :]
		requiredAAs = requiredAAs[inversePermutationIndexes, :]
		updatedAAs = updatedAAs[inversePermutationIndexes, :]

		assert np.all(updatedAAs <= requiredAAs), "Polypeptides got elongated more than possible!"

		updatedMass = massDiffProtein + np.dot(
			(updatedAAs - assignedAAs), self.aaWeightsIncorporated
			).flatten()


		didInitialize = (
			(assignedAAs.sum(axis = 1) == 0) &
			(updatedAAs.sum(axis = 1) > 0)
			)

		activeRibosomes.attrIs(
			assignedAAs = updatedAAs,
			massDiffProtein = updatedMass
			)

		terminatedProteins = np.zeros_like(self.bulkMonomers.counts())

		didTerminate = (requiredAAs == updatedAAs).all(axis = 1)

		for moleculeIndex, molecule in enumerate(activeRibosomes):
			if didTerminate[moleculeIndex]:
				terminatedProteins[molecule.attr('proteinIndex')] += 1
				self.activeRibosomes.moleculeDel(molecule)

		nTerminated = didTerminate.sum()
		nInitialized = didInitialize.sum()
		nElongations = aasUsed.sum()

		self.aas.countsDec(aasUsed)

		self.bulkMonomers.countsIs(terminatedProteins)

		self.ribosomeSubunits.countsInc(nTerminated)

		self.h2o.countInc(nElongations)


def getWorkAssignment(dataSize, thisTask, totalTasks):
	startPos = sum(dataSize // totalTasks + (i < (dataSize % totalTasks)) for i in xrange(thisTask))
	nElements = dataSize // totalTasks + (thisTask < (dataSize % totalTasks))

	return startPos, nElements