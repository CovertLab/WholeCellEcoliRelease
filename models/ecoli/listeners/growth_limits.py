"""
GrowthLimits
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

# from numpy.lib.recfunctions import merge_arrays

VERBOSE = False

class GrowthLimits(wholecell.listeners.listener.Listener):
	""" GrowthLimits """

	_name = 'GrowthLimits'

	# Constructor
	def __init__(self, *args, **kwargs):
		super(GrowthLimits, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, sim_data):
		super(GrowthLimits, self).initialize(sim, sim_data)

		# Computed, saved attributes
		self.aaIds = sim_data.molecule_groups.amino_acids
		self.ntpIds = sim_data.molecule_groups.ntps
		self.uncharged_trna_ids = sim_data.process.transcription.rna_data['id'][sim_data.process.transcription.rna_data['is_tRNA']].tolist()
		self.charged_trna_ids = sim_data.process.transcription.charged_trna_names
		self.aa_enzymes = list(sim_data.process.metabolism.aa_enzymes)
		self.aa_exporter_names = list(sim_data.process.metabolism.aa_export_transporters_names)

	# Allocate memory
	def allocate(self):
		super(GrowthLimits, self).allocate()

		n_aa = len(self.aaIds)
		n_enzymes = len(self.aa_enzymes)
		n_exporters = len(self.aa_exporter_names)

		# For translation
		self.activeRibosomeAllocated = 0

		self.aaPoolSize = np.zeros(n_aa, np.float64)
		self.aaRequestSize = np.zeros(n_aa, np.float64)
		self.aaAllocated = np.zeros(n_aa, np.float64)
		self.aasUsed = np.zeros(n_aa, np.float64)

		# For charging function
		self.synthetase_conc = np.zeros(n_aa, np.float64)
		self.uncharged_trna_conc = np.zeros(n_aa, np.float64)
		self.charged_trna_conc = np.zeros(n_aa, np.float64)
		self.aa_conc = np.zeros(n_aa, np.float64)
		self.ribosome_conc = 0.
		self.fraction_aa_to_elongate = np.zeros(n_aa, np.float64)

		# Charging results
		n_uncharged_trna = len(self.uncharged_trna_ids)
		self.fraction_trna_charged = np.zeros(n_uncharged_trna, np.float64)
		self.net_charged = np.zeros(n_uncharged_trna, int)

		# For transcription
		n_ntp = len(self.ntpIds)
		self.ntpPoolSize = np.zeros(n_ntp, np.float64)
		self.ntpRequestSize = np.zeros(n_ntp, np.float64)
		self.ntpAllocated = np.zeros(n_ntp, np.float64)
		self.ntpUsed = np.zeros(n_ntp, np.float64)

		self.ppgpp_conc = 0.
		self.rela_conc = 0.
		self.spot_conc = 0.
		self.rela_syn = 0.
		self.spot_syn = 0.
		self.spot_deg = 0.

		self.original_aa_supply = np.zeros(n_aa, np.float64)
		self.aa_supply = np.zeros(n_aa, np.float64)
		self.aa_synthesis = np.zeros(n_aa, np.float64)
		self.aa_import = np.zeros(n_aa, np.float64)
		self.aa_export = np.zeros(n_aa, np.float64)
		self.aa_supply_enzymes = np.zeros(n_enzymes, int)
		self.aa_exporters = np.zeros(n_exporters, int)
		self.aa_supply_aa_conc = np.zeros(n_aa, np.float64)
		self.aa_supply_fraction = np.zeros(n_aa, np.float64)
		self.aa_in_media = np.zeros(n_aa, bool)

		self.aaCountDiff = np.zeros(n_aa, np.float64)
		self.trnaCharged = np.zeros(n_aa, np.float64)

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'aaPoolSize': 'aaIds',
			'aaRequestSize': 'aaIds',
			'aaAllocated': 'aaIds',
			'aasUsed': 'aaIds',
			'synthetase_conc': 'aaIds',
			'uncharged_trna_conc': 'aaIds',
			'charged_trna_conc': 'aaIds',
			'aa_conc': 'aaIds',
			'fraction_aa_to_elongate': 'aaIds',
			'fraction_trna_charged': 'uncharged_trna_ids',
			'net_charged': 'uncharged_trna_ids',
			'ntpPoolSize': 'ntpIds',
			'ntpRequestSize': 'ntpIds',
			'ntpAllocated': 'ntpIds',
			'ntpUsed': 'ntpIds',
			'original_aa_supply': 'aaIds',
			'aa_supply': 'aaIds',
			'aa_synthesis': 'aaIds',
			'aa_import': 'aaIds',
			'aa_export': 'aaIds',
			'aa_supply_enzymes': 'aa_enzymes',
			'aa_exporters': 'aa_exporter_names',
			'aa_supply_aa_conc': 'aaIds',
			'aa_supply_fraction': 'aaIds',
			'aa_in_media': 'aaIds',
			'trnaCharged':'aaIds',
			'aaCountDiff':'aaIds'
			}

		tableWriter.writeAttributes(
			aaIds = self.aaIds,
			uncharged_trna_ids = self.uncharged_trna_ids,
			ntpIds = self.ntpIds,
			subcolumns = subcolumns,
			aa_enzymes = self.aa_enzymes,
			aa_exporter_names = self.aa_exporter_names,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			activeRibosomeAllocated = self.activeRibosomeAllocated,
			aaPoolSize = self.aaPoolSize,
			aaRequestSize = self.aaRequestSize,
			aaAllocated = self.aaAllocated,
			aasUsed = self.aasUsed,
			synthetase_conc = self.synthetase_conc,
			uncharged_trna_conc = self.uncharged_trna_conc,
			charged_trna_conc = self.charged_trna_conc,
			aa_conc = self.aa_conc,
			ribosome_conc = self.ribosome_conc,
			fraction_aa_to_elongate = self.fraction_aa_to_elongate,
			fraction_trna_charged = self.fraction_trna_charged,
			net_charged = self.net_charged,
			ntpPoolSize = self.ntpPoolSize,
			ntpRequestSize = self.ntpRequestSize,
			ntpAllocated = self.ntpAllocated,
			ntpUsed = self.ntpUsed,
			ppgpp_conc = self.ppgpp_conc,
			rela_conc = self.rela_conc,
			spot_conc = self.spot_conc,
			rela_syn = self.rela_syn,
			spot_syn = self.spot_syn,
			spot_deg = self.spot_deg,
			original_aa_supply = self.aa_supply,
			aa_supply = self.aa_supply,
			aa_synthesis = self.aa_synthesis,
			aa_import = self.aa_import,
			aa_export = self.aa_export,
			aa_supply_enzymes = self.aa_supply_enzymes,
			aa_exporters = self.aa_exporters,
			aa_supply_aa_conc = self.aa_supply_aa_conc,
			aa_supply_fraction = self.aa_supply_fraction,
			aa_in_media = self.aa_in_media,
			aaCountDiff = self.aaCountDiff,
			trnaCharged = self.trnaCharged,
			)
