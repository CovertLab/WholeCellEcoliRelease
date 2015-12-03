#!/usr/bin/env python

"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- complexes

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/16/2013
"""

from __future__ import division

import numpy as np

import wholecell.processes.process
from wholecell.utils.constants import REQUEST_PRIORITY_DEGRADATION
from wholecell.utils import units

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
	def initialize(self, sim, sim_data):
		super(ProteinDegradation, self).initialize(sim, sim_data)

		# Metabolite IDs for S matrix

		h2oId = ["WATER[c]"]

		metaboliteIds = sim_data.moleculeGroups.aaIDs + h2oId

		aaIdxs = np.arange(0, len(sim_data.moleculeGroups.aaIDs))
		h2oIdx = metaboliteIds.index('WATER[c]')

		# Protein IDs for S matrix
		proteinIds = sim_data.process.translation.monomerData['id']

		# Proteins
		self.proteinDegRates = sim_data.process.translation.monomerData['degRate'].asNumber(1 / units.s) * self.timeStepSec

		self.proteinLengths = sim_data.process.translation.monomerData['length']

		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIdxs, :] = np.transpose(sim_data.process.translation.monomerData["aaCounts"].asNumber())
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIdxs, :], axis = 0) - 1)

		# Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView('WATER[c]')
		self.proteins = self.bulkMoleculesView(proteinIds)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

		# TODO: Curate which proteases to use here
		# self.protease = self.bulkMoleculeView('EG11259-MONOMER[c]')

	def calculateRequest(self):
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self.proteinDegRates * self.proteins.total()),
			self.proteins.total()
			)

		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

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
