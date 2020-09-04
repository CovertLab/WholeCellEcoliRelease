"""
SimulationData for replication process

@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 02/13/2015
"""

from __future__ import absolute_import, division, print_function

import numpy as np
import collections

from wholecell.sim.simulation import MAX_TIME_STEP
from wholecell.utils import units
from wholecell.utils.polymerize import polymerize
from wholecell.utils.random import stochasticRound


PROCESS_MAX_TIME_STEP = 2.


class Replication(object):
	"""
	SimulationData for the replication process
	"""

	def __init__(self, raw_data, sim_data):
		self.max_time_step = min(MAX_TIME_STEP, PROCESS_MAX_TIME_STEP)

		self._n_nt_types = len(sim_data.dntp_code_to_id_ordered)
		self._c_period = sim_data.growthRateParameters.c_period.asNumber(units.min)
		self._d_period = sim_data.growthRateParameters.d_period.asNumber(units.min)

		self._buildSequence(raw_data, sim_data)
		self._buildGeneData(raw_data, sim_data)
		self._buildReplication(raw_data, sim_data)
		self._buildMotifs(raw_data, sim_data)
		self._build_elongation_rates(raw_data, sim_data)

	def _buildSequence(self, raw_data, sim_data):
		self.genome_sequence = raw_data.genome_sequence
		self.genome_sequence_rc = self.genome_sequence.reverse_complement()
		self.genome_length = len(self.genome_sequence)
		self.genome_A_count = self.genome_sequence.count("A")
		self.genome_T_count = self.genome_sequence.count("T")
		self.genome_G_count = self.genome_sequence.count("G")
		self.genome_C_count = self.genome_sequence.count("C")

	def _buildGeneData(self, raw_data, sim_data):
		"""
		Build gene-associated simulation data from raw data.
		"""

		def extract_data(raw, key):
			data = [row[key] for row in raw]
			dtype = 'U{}'.format(max(len(d) for d in data if d is not None))
			return data, dtype

		names, name_dtype = extract_data(raw_data.genes, 'id')
		symbols, symbol_dtype = extract_data(raw_data.genes, 'symbol')
		rna_ids, rna_dtype = extract_data(raw_data.genes, 'rnaId')

		self.geneData = np.zeros(
			len(raw_data.genes),
			dtype=[('name', name_dtype),
				('symbol', symbol_dtype),
				('rnaId', rna_dtype)])

		self.geneData['name'] = names
		self.geneData['symbol'] = symbols
		self.geneData['rnaId'] = rna_ids

	def _buildReplication(self, raw_data, sim_data):
		"""
		Build replication-associated simulation data from raw data.
		"""
		# Map ATGC to 8 bit integers
		numerical_sequence = np.empty(self.genome_length, np.int8)
		ntMapping = collections.OrderedDict(
			[(ntp_code, i) for i, ntp_code in enumerate(sim_data.dntp_code_to_id_ordered.keys())])
		for i, letter in enumerate(raw_data.genome_sequence):
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

		self.forward_complement_sequence = self._get_complement_sequence(self.forward_sequence)
		self.reverse_complement_sequence = self._get_complement_sequence(self.reverse_sequence)

		assert self.forward_sequence.size + self.reverse_sequence.size == self.genome_length

		# Log array of lengths of each replichore
		self.replichore_lengths = np.array([
			self.forward_sequence.size,
			self.reverse_sequence.size,
			])

		# Determine size of the matrix used by polymerize function
		maxLen = np.int64(
			self.replichore_lengths.max()
			+ self.max_time_step * sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(units.nt / units.s)
		)

		self.replication_sequences = np.empty((4, maxLen), np.int8)
		self.replication_sequences.fill(polymerize.PAD_VALUE)

		self.replication_sequences[0, :self.forward_sequence.size] = self.forward_sequence
		self.replication_sequences[1, :self.forward_complement_sequence.size] = self.forward_complement_sequence
		self.replication_sequences[2, :self.reverse_sequence.size] = self.reverse_sequence
		self.replication_sequences[3, :self.reverse_complement_sequence.size] = self.reverse_complement_sequence

		# Get polymerized nucleotide weights
		self.replicationMonomerWeights = (
			(sim_data.getter.getMass(sim_data.moleculeGroups.dntps)
			- sim_data.getter.getMass([sim_data.moleculeIds.ppi]))
			/ sim_data.constants.nAvogadro
		)

		# Placeholder value for "child_domains" attribute of domains without
		# children domains
		self.no_child_place_holder = -1

	def _buildMotifs(self, raw_data, sim_data):
		"""
		Build simulation data associated with sequence motifs from raw_data.
		Coordinates of all motifs are calculated based on the given sequences
		of the genome and the motifs.
		"""
		# Get coordinates of oriC's and terC's
		self.oric_coordinates = raw_data.parameters['oriCCenter'].asNumber()
		self.terc_coordinates = raw_data.parameters['terCCenter'].asNumber()

		# Initialize dictionary of motif coordinates
		self.motif_coordinates = dict()

		for motif in raw_data.sequence_motifs:
			# Keys are the IDs of the motifs; Values are arrays of the motif's
			# coordinates.
			self.motif_coordinates[motif["id"]] = self._get_motif_coordinates(
				motif["length"], motif["sequences"])

	def _get_complement_sequence(self, sequenceVector):
		"""
		Calculates the vector for a complement sequence of a DNA sequence given
		in vector form.
		"""
		return (self._n_nt_types - 1) - sequenceVector

	def _get_motif_coordinates(self, motif_length, motif_sequences):
		"""
		Finds the coordinates of all sequence motifs of a specific type. The
		coordinates are given as the positions of the midpoint of the motif
		relative to the oriC in base pairs.
		"""
		# Append the first n bases to the end of the sequences to account for
		# motifs that span the two endpoints
		extended_sequence = self.genome_sequence + self.genome_sequence[:(motif_length - 1)]
		extended_rc_sequence = self.genome_sequence_rc + self.genome_sequence_rc[:(motif_length - 1)]

		# Initialize list for coordinates of motifs
		coordinates = []

		# Loop through all possible motif sequences
		for sequence in motif_sequences:
			# Find occurrences of the motif in the original sequence
			loc = extended_sequence.find(sequence) # Returns -1 if not found

			while loc != -1:
				coordinates.append(loc + motif_length//2)
				loc = extended_sequence.find(sequence, loc + 1)

			# Find occurrences of the motif in the reverse complement
			loc = extended_rc_sequence.find(sequence)

			while loc != -1:
				coordinates.append(
					self.genome_length - loc - motif_length + motif_length//2)
				loc = extended_rc_sequence.find(sequence, loc + 1)

		# Compute coordinates relative to oriC
		motif_coordinates = self._get_relative_coordinates(np.array(coordinates))
		motif_coordinates.sort()

		return motif_coordinates

	def _get_relative_coordinates(self, coordinates):
		"""
		Converts an array of genomic coordinates into coordinates relative to
		the origin of replication.
		"""
		relative_coordinates = (
			(coordinates - self.terc_coordinates)
			% self.genome_length + self.terc_coordinates - self.oric_coordinates)

		relative_coordinates[relative_coordinates < 0] += 1

		return relative_coordinates

	def _build_elongation_rates(self, raw_data, sim_data):
		self.basal_elongation_rate = int(
			round(sim_data.growthRateParameters.dnaPolymeraseElongationRate.asNumber(
			units.nt / units.s)))

	def make_elongation_rates(self, random, replisomes, base, time_step):
		rates = np.full(
			replisomes,
			stochasticRound(random, base * time_step),
			dtype=np.int64)

		return rates

	def get_average_copy_number(self, tau, coords):
		"""
		Calculates the average copy number of a gene throughout the cell cycle
		given the location of the gene in coordinates.

		Args:
			tau (float): expected doubling time in minutes
			coords (int or ndarray[int]): chromosome coordinates of genes

		Returns:
			float or ndarray[float] (matches length of coords): average copy
			number of each gene expected at a doubling time, tau
		"""

		right_replichore_length = self.replichore_lengths[0]
		left_replichore_length = self.replichore_lengths[1]

		# Calculate the relative position of the gene along the chromosome
		# from its coordinate
		relative_pos = np.array(coords, float)
		relative_pos[coords > 0] = relative_pos[coords > 0] / right_replichore_length
		relative_pos[coords < 0] = -relative_pos[coords < 0] / left_replichore_length

		# Return the predicted average copy number
		n_avg_copy = 2**(((1 - relative_pos) * self._c_period + self._d_period) / tau)
		return n_avg_copy
