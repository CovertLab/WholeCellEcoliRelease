"""
SimulationData for replication process

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np
import collections

from wholecell.utils import units
from wholecell.utils.polymerize import polymerize

MAX_TIMESTEP_LEN = 2

class Replication(object):
	"""
	SimulationData for the replication process
	"""

	def __init__(self, raw_data, sim_data):
		self._n_nt_types = len(sim_data.dNtpOrder)

		self._buildSequence(raw_data, sim_data)
		self._buildGeneData(raw_data, sim_data)
		self._buildReplication(raw_data, sim_data)

	def _buildSequence(self, raw_data, sim_data):
		self.genome_sequence = raw_data.genome_sequence
		self.genome_length = len(self.genome_sequence)
		self.genome_A_count = self.genome_sequence.count("A")
		self.genome_T_count = self.genome_sequence.count("T")
		self.genome_G_count = self.genome_sequence.count("G")
		self.genome_C_count = self.genome_sequence.count("C")

	def _buildGeneData(self, raw_data, sim_data):
		"""
		Build gene-associated simulation data from raw data.
		"""

		self.geneData = np.zeros(
			len(raw_data.genes),
			dtype=[('name', 'a50'),
				('symbol', 'a7'),
				('rnaId', 'a50'),
				('monomerId', 'a50')])

		self.geneData['name'] = [x['id'] for x in raw_data.genes]
		self.geneData['symbol'] = [x['symbol'] for x in raw_data.genes]
		self.geneData['rnaId'] = [x['rnaId'] for x in raw_data.genes]
		self.geneData['monomerId'] = [x['monomerId'] for x in raw_data.genes]

	def _buildReplication(self, raw_data, sim_data):
		"""
		Build replication-associated simulation data from raw data.
		"""
		# Map ATGC to 8 bit integers
		numerical_sequence = np.empty(self.genome_length, np.int8)
		ntMapping = collections.OrderedDict([(ntpId, i) for i, ntpId in enumerate(sim_data.dNtpOrder)])
		for i,letter in enumerate(raw_data.genome_sequence):
			numerical_sequence[i] = ntMapping[letter] # Build genome sequence as small integers

		# Create 4 possible polymerization sequences
		oric_coordinate = raw_data.parameters['oriCCenter'].asNumber()
		terc_coordinate = raw_data.parameters['terCCenter'].asNumber()

		# Forward sequence includes oriC
		self.forward_sequence = numerical_sequence[
			np.hstack((np.arange(oric_coordinate, self.genome_length),
			np.arange(0, terc_coordinate)))]

		# Reverse sequence includes terC
		self.reverse_sequence = numerical_sequence[
			np.arange(oric_coordinate - 1, terc_coordinate - 1, -1)]

		self.forward_complement_sequence = self._reverseComplement(self.forward_sequence)
		self.reverse_complement_sequence = self._reverseComplement(self.reverse_sequence)

		assert self.forward_sequence.size + self.reverse_sequence.size == self.genome_length

		# Log array of lengths of each replichore
		self.replichore_lengths = np.array([
			self.forward_sequence.size,
			self.reverse_sequence.size,
			])

		# Determine size of the matrix used by polymerize function
		maxLen = np.int64(
			self.replichore_lengths.max()
			+ MAX_TIMESTEP_LEN * sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(units.nt / units.s)
		)

		self.replication_sequences = np.empty((4, maxLen), np.int8)
		self.replication_sequences.fill(polymerize.PAD_VALUE)

		self.replication_sequences[0, :self.forward_sequence.size] = self.forward_sequence
		self.replication_sequences[1, :self.forward_complement_sequence.size] = self.forward_complement_sequence
		self.replication_sequences[2, :self.reverse_sequence.size] = self.reverse_sequence
		self.replication_sequences[3, :self.reverse_complement_sequence.size] = self.reverse_complement_sequence

		# Get polymerized nucleotide weights
		self.replicationMonomerWeights = (
			(sim_data.getter.getMass(sim_data.moleculeGroups.dNtpIds)
			- sim_data.getter.getMass(["PPI[c]"]))
			/ raw_data.constants['nAvogadro']
		)

	def _reverseComplement(self, sequenceVector):
		return (self._n_nt_types - 1) - sequenceVector
