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
		# Constants
		self.proteinLengths = None		# Protein lengths
		self.proteinDegRates = None		# Protein degradation rates (1/s)
		self.proteinDegSMatrix = None	# Protein degradation stoichiometry matrix [metabolite x rna]

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, kb):
		super(ProteinDegradation, self).initialize(sim, kb)

		# Metabolite IDs for S matrix
		aaIds = [
		"ALA-L[c]", "ARG-L[c]", "ASN-L[c]", "ASP-L[c]", "CYS-L[c]", "GLU-L[c]", "GLN-L[c]", "GLY[c]",
		"HIS-L[c]", "ILE-L[c]", "LEU-L[c]", "LYS-L[c]", "MET-L[c]", "PHE-L[c]", "PRO-L[c]", "SELNP[c]",
		 "SER-L[c]", "THR-L[c]", "TRP-L[c]", "TYR-L[c]", "VAL-L[c]"
		]
		h2oId = ["H2O[c]"]

		metaboliteIds = aaIds + h2oId

		self._aaIds = np.arange(0, 21)
		self._h2oIdx = metaboliteIds.index('H2O[c]')

		# Protein IDs for S matrix
		proteinIds = kb.monomerData['id']

		# Proteins
		# TODO: Generate degradation rates
		self.proteinDegRates = kb.monomerData['degRate']

		self.proteinLengths = kb.monomerData['length']

		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[self._aaIds, :] = np.transpose(kb.monomerData['aaCounts'])
		self.proteinDegSMatrix[self._h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[self._aaIds, :], axis = 0) - 1)

		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView('H2O[c]')
		self.proteins = self.bulkMoleculesView(proteinIds)

		# TODO: Curate which proteases to use here
		# self.protease = self.bulkMoleculeView('EG11259-MONOMER[c]')

	def calculateRequest(self):
		nProteinsToDegrade = np.fmin(
			self.randStream.poissrnd(self.proteinDegRates * self.proteins.total() * self.timeStepSec),
			self.proteins.total()
			)

		nReactions = np.dot(self.proteinLengths, nProteinsToDegrade)

		# TODO: length-1 reactions? for water not nReactions. Check rna degradation as well.
		self.h2o.requestIs(nReactions)
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
		self.h2o.countsIs(0)
