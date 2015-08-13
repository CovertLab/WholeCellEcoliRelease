"""
SimulationData for replication process

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import division

import numpy as np
import collections

from wholecell.utils import units

class Replication(object):
	""" Replication """

	def __init__(self, raw_data, sim_data):
		self._n_nt_types = len(sim_data.dNtpOrder)

		self._buildSequence(raw_data, sim_data)
		self._buildGeneData(raw_data, sim_data)
		self._buildReplication(raw_data, sim_data)
		self._determineIniationCellMasses(raw_data, sim_data)

	def _buildSequence(self, raw_data, sim_data):
		self.genome_sequence = raw_data.genome_sequence
		self.genome_length = len(self.genome_sequence)
		self.genome_A_count = self.genome_sequence.count("A")
		self.genome_T_count = self.genome_sequence.count("T")
		self.genome_G_count = self.genome_sequence.count("G")
		self.genome_C_count = self.genome_sequence.count("C")

	def _buildGeneData(self, raw_data, sim_data):
		genomeLength = len(raw_data.genome_sequence)
		self.geneData = np.zeros(len(raw_data.genes),
			dtype = [('name'				,	'a50'),
					('rnaId'                ,   'a50'),
					('endCoordinate'		,	'int64')])

		self.geneData['name'] = [x['id'] for x in raw_data.genes]
		self.geneData['rnaId'] = [x['rnaId'] for x in raw_data.genes]
		self.geneData['endCoordinate'] = [(x['coordinate'] + x['length']) % genomeLength if x['direction'] == '+' else (x['coordinate'] - x['length']) % genomeLength for x in raw_data.genes]

	def _buildReplication(self, raw_data, sim_data):
		# Map ATGC to 8 bit integers
		numerical_sequence = np.empty(self.genome_length, np.int8)
		ntMapping = collections.OrderedDict([(ntpId, i) for i, ntpId in enumerate(sim_data.dNtpOrder)])
		for i,letter in enumerate(raw_data.genome_sequence):
			numerical_sequence[i] = ntMapping[letter] # Build genome sequence as small integers

		# Create 4 possible polymerization sequences
		oriC = raw_data.parameters['oriCCenter'].asNumber()
		terC = raw_data.parameters['terCCenter'].asNumber()

		self.forward_sequence = numerical_sequence[np.hstack((np.arange(oriC, self.genome_length),np.arange(0, terC)))]
		self.reverse_sequence = numerical_sequence[np.arange(oriC, terC, -1)]
		self.forward_complement_sequence = self._reverseComplement(self.forward_sequence)
		self.reverse_complement_sequence = self._reverseComplement(self.reverse_sequence)

		assert self.forward_sequence.size + self.reverse_sequence.size == self.genome_length

		# Build sequence matrix for polymerize function in replication process
		self.sequence_lengths = np.array([
										self.forward_sequence.size,
										self.reverse_sequence.size,
										self.forward_complement_sequence.size,
										self.reverse_complement_sequence.size
										])

		maxLen = np.int64(
			self.sequence_lengths.max()
			+ sim_data.constants.dnaPolymeraseElongationRate.asNumber(units.nt / units.s)# * sim_data.timeStepSec # TODO: FIX
			)

		from wholecell.utils.polymerize import PAD_VALUE
		self.replication_sequences = np.empty((4, maxLen), np.int8)
		self.replication_sequences.fill(PAD_VALUE)

		self.replication_sequences[0, :self.forward_sequence.size] = self.forward_sequence
		self.replication_sequences[1, :self.reverse_sequence.size] = self.reverse_sequence
		self.replication_sequences[2, :self.forward_complement_sequence.size] = self.forward_complement_sequence
		self.replication_sequences[3, :self.reverse_complement_sequence.size] = self.reverse_complement_sequence

		# Get polymeriezd nucleotide weights
		self.replicationMonomerWeights = (
			(
				sim_data.getter.getMass(sim_data.moleculeGroups.polymerizedDNT_IDs)
				- sim_data.getter.getMass(["PPI[c]"])
				)
			/ raw_data.constants['nAvogadro']
			)

	def _reverseComplement(self, sequenceVector):
		return (self._n_nt_types - 1) - sequenceVector


	def _determineIniationCellMasses(self, raw_data, sim_data):
		if sim_data.doubling_time >= 60. * units.min:
			self.cellMassReplicationInitiation = 1 * sim_data.mass.avgCell60MinDoublingTimeTotalMassInit
		elif sim_data.doubling_time >= 30. * units.min:
			self.cellMassReplicationInitiation = 2 * sim_data.mass.avgCell60MinDoublingTimeTotalMassInit
		elif sim_data.doubling_time >= 20. * units.min:
			self.cellMassReplicationInitiation = 4 * sim_data.mass.avgCell60MinDoublingTimeTotalMassInit
		else:
			raise Exception("Un-expected doubling time!")
