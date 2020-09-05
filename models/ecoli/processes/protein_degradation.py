"""
ProteinDegradation

Protein degradation sub-model. Encodes molecular simulation of protein degradation as a Poisson process

TODO:
- protein complexes
- add protease functionality

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 5/16/2013
"""

from __future__ import absolute_import, division, print_function

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
		self.proteinDegSMatrix = None	# Protein degradation stoichiometry matrix [metabolite x rna]

		# Views
		self.metabolites = None
		self.h2o = None
		self.proteins = None

		super(ProteinDegradation, self).__init__()

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(ProteinDegradation, self).initialize(sim, sim_data)

		# Load protein degradation rates (based on N-end rule)
		self.rawDegRate = sim_data.process.translation.monomer_data['deg_rate'].asNumber(1 / units.s)

		shuffleIdxs = None
		if hasattr(sim_data.process.translation, "monomerDegRateShuffleIdxs") and sim_data.process.translation.monomerDegRateShuffleIdxs is not None:
			shuffleIdxs = sim_data.process.translation.monomerDegRateShuffleIdxs
			self.rawDegRate = self.rawDegRate[shuffleIdxs]

		# Build metabolite IDs for S matrix
		h2oId = [sim_data.molecule_ids.water]
		metaboliteIds = sim_data.molecule_groups.amino_acids + h2oId
		aaIdxs = np.arange(0, len(sim_data.molecule_groups.amino_acids))
		h2oIdx = metaboliteIds.index(sim_data.molecule_ids.water)

		# Build protein IDs for S matrix
		proteinIds = sim_data.process.translation.monomer_data['id']

		# Load protein length
		self.proteinLengths = sim_data.process.translation.monomer_data['length']

		# Build S matrix
		self.proteinDegSMatrix = np.zeros((len(metaboliteIds), len(proteinIds)), np.int64)
		self.proteinDegSMatrix[aaIdxs, :] = np.transpose(sim_data.process.translation.monomer_data['aa_counts'].asNumber())
		self.proteinDegSMatrix[h2oIdx, :]  = -(np.sum(self.proteinDegSMatrix[aaIdxs, :], axis = 0) - 1)

		# Build Views
		self.metabolites = self.bulkMoleculesView(metaboliteIds)
		self.h2o = self.bulkMoleculeView(sim_data.molecule_ids.water)
		self.proteins = self.bulkMoleculesView(proteinIds)

		self.bulkMoleculesRequestPriorityIs(REQUEST_PRIORITY_DEGRADATION)

	def calculateRequest(self):

		# Determine how many proteins to degrade based on the degradation rates and counts of each protein
		nProteinsToDegrade = np.fmin(
			self.randomState.poisson(self._proteinDegRates() * self.proteins.total_counts()),
			self.proteins.total_counts()
			)

		# Determine the number of hydrolysis reactions
		nReactions = np.dot(self.proteinLengths.asNumber(), nProteinsToDegrade)

		# Determine the amount of water required to degrade the selected proteins
		# Assuming one N-1 H2O is required per peptide chain length N
		self.h2o.requestIs(nReactions - np.sum(nProteinsToDegrade))
		self.proteins.requestIs(nProteinsToDegrade)


	def evolveState(self):

		# Degrade selected proteins, release amino acids from those proteins back into the cell, 
		# and consume H_2O that is required for the degradation process
		self.metabolites.countsInc(np.dot(
			self.proteinDegSMatrix,
			self.proteins.counts()
			))
		self.proteins.countsIs(0)

	def _proteinDegRates(self):
		return self.rawDegRate * self.timeStepSec()
