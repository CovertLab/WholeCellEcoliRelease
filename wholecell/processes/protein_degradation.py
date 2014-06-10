#!/usr/bin/env python

"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/16/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process


class ProteinDegradation(wholecell.processes.process.Process):
	""" ProteinDegradation """

	_name = "ProteinDegradation"

	# Constructor
	def __init__(self):
		# Parameters
		self.proteinLengths = None		# Protein lengths
		self.proteinDegRates = None		# Protein degradation rates (1/s)
		self.proteinDegSMatrix = None	# Protein degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.h2o = None
		self.proteins = None

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinDegradation, self).initialize(sim, kb)

		# Metabolite IDs for S matrix
		aaIds = kb.aaIDs[:]

		# TODO: Remove hack of deleting selenocysteine this way
		selenocysteineIdx = aaIds.index("SEC-L[c]")
		del aaIds[selenocysteineIdx]

		proteinAACounts = np.delete(
			kb.monomerData["aaCounts"].magnitude, selenocysteineIdx, 1
			)

		h2oId = ["H2O[c]"]

		metaboliteIds = aaIds + h2oId

		aaIds = np.arange(0, len(aaIds))
		h2oIdx = metaboliteIds.index('H2O[c]')

		# Protein IDs for S matrix
		proteinIds = kb.monomerData['id']

		# Proteins
		self.proteinDegRates = kb.monomerData['degRate'].to('1 / s').magnitude

		self.proteinLengths = kb.monomerData['length']

		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIds, :] = np.transpose(proteinAACounts)
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIds, :], axis = 0) - 1)

		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proteins = self.bulkMoleculesView(proteinIds)

		# TODO: Curate which proteases to use here
		# self.protease = self.bulkMoleculeView('EG11259-MONOMER[c]')

	def calculateRequest(self):
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self.proteinDegRates * self.proteins.total() * self.timeStepSec),
			self.proteins.total()
			)

		nReactions = np.dot(self.proteinLengths, nProteinsToDegrade)

		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)
		#self.protease.requestAll()
		

	def evolveState(self):
		# Check if protease expressed
		#if self.protease.count() == 0:
		#	return

		# Degrade protein
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))

		self.proteins.countsIs(0)
		self.h2o.countIs(0)
