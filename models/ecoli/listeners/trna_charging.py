"""
TrnaCharging
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener
from wholecell.utils import units

class TrnaCharging(wholecell.listeners.listener.Listener):
	""" TrnaCharging """

	_name = 'TrnaCharging'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(TrnaCharging, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TrnaCharging, self).initialize(sim, sim_data)

		rna_data = sim_data.process.transcription.rna_data
		trnas = rna_data['id'][rna_data['is_tRNA']]
		self.n_trnas = len(trnas)
		self.n_codons = len(sim_data.relation.codons)
		self.trna_codon_pairs = sim_data.relation.trna_codon_pairs
		self.n_trna_codon_pairs = len(self.trna_codon_pairs)
		self.n_amino_acids = len(sim_data.molecule_groups.amino_acids)


		# Arginine biosynthesis
		monomer_data = sim_data.process.translation.monomer_data
		argA_index = np.where(monomer_data['id'] == 'N-ACETYLTRANSFER-MONOMER[c]')[0][0]

		# add 1 to represent terminated ribosomes
		self.argA_listener_size = int(1 + monomer_data['length'][argA_index].asNumber(units.aa))

	# Allocate memory
	def allocate(self):
		super(TrnaCharging, self).allocate()

		self.charging_events = np.zeros(self.n_trnas, np.int64)
		self.cleaved = 0
		self.codons_read = np.zeros(self.n_codons, np.int64)
		self.codons_to_trnas_counter = np.zeros(self.n_trna_codon_pairs, np.int64)
		self.collision_removed_ribosomes_argA = np.zeros(self.argA_listener_size, np.int64)
		self.initial_disagreements = np.zeros(self.n_codons)
		self.initiated = 0
		self.not_cleaved = 0
		self.reading_events = np.zeros(self.n_trnas, np.int64)
		self.ribosome_positions_argA = np.zeros(self.argA_listener_size, np.int64)
		self.ribosome_initiation_argA = np.array([], np.int64)

		self.saturation_trna = np.zeros(self.n_trnas)
		self.turnover = np.zeros(self.n_amino_acids)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			trna_codon_pairs = self.trna_codon_pairs,
			)
		tableWriter.set_variable_length_columns(
			'ribosome_initiation_argA',
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			charging_events = self.charging_events,
			cleaved = self.cleaved,
			codons_read = self.codons_read,
			codons_to_trnas_counter = self.codons_to_trnas_counter,
			collision_removed_ribosomes_argA = self.collision_removed_ribosomes_argA,
			initial_disagreements = self.initial_disagreements,
			initiated = self.initiated,
			not_cleaved = self.not_cleaved,
			reading_events = self.reading_events,
			ribosome_positions_argA = self.ribosome_positions_argA,
			ribosome_initiation_argA = self.ribosome_initiation_argA,
			saturation_trna = self.saturation_trna,
			turnover = self.turnover,
			)
