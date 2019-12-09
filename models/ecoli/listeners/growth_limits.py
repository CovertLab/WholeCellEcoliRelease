#!/usr/bin/env python

"""
GrowthLimits

@author: Nick Ruggero
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 6/25/15
"""

from __future__ import division

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
		self.aaIds = sim_data.moleculeGroups.aaIDs
		self.ntpIds = sim_data.moleculeGroups.ntpIds
		self.uncharged_trna_ids = sim_data.process.transcription.rnaData['id'][sim_data.process.transcription.rnaData['isTRna']].tolist()
		self.charged_trna_ids = sim_data.process.transcription.charged_trna_names
		self.n_aas = len(self.aaIds)

	# Allocate memory
	def allocate(self):
		super(GrowthLimits, self).allocate()

		# For translation
		self.gtpPoolSize = 0
		self.gtpRequestSize = 0
		self.gtpAllocated = 0
		self.gtpUsed = 0

		self.activeRibosomeAllocated = 0

		n_aa = len(self.aaIds)
		self.aaPoolSize = np.zeros(n_aa, np.float64)
		self.aaRequestSize = np.zeros(n_aa, np.float64)
		self.aaAllocated = np.zeros(n_aa, np.float64)
		self.aasUsed = np.zeros(n_aa, np.float64)

		n_uncharged_trna = len(self.uncharged_trna_ids)
		self.fraction_trna_charged = np.zeros(n_uncharged_trna, np.float64)
		self.net_charged = np.zeros(n_uncharged_trna, np.int)

		# For transcription
		n_ntp = len(self.ntpIds)
		self.ntpPoolSize = np.zeros(n_ntp, np.float64)
		self.ntpRequestSize = np.zeros(n_ntp, np.float64)
		self.ntpAllocated = np.zeros(n_ntp, np.float64)
		self.ntpUsed = np.zeros(n_ntp, np.float64)

		self.rela_syn = 0.
		self.spot_syn = 0.
		self.spot_deg = 0.

	def update(self):
		pass

	def tableCreate(self, tableWriter):
		subcolumns = {
			'aaPoolSize': 'aaIds',
			'aaRequestSize': 'aaIds',
			'aaAllocated': 'aaIds',
			'aasUsed': 'aaIds',

			'fraction_trna_charged': 'uncharged_trna_ids',
			'net_charged': 'uncharged_trna_ids',

			'ntpPoolSize': 'ntpIds',
			'ntpRequestSize': 'ntpIds',
			'ntpAllocated': 'ntpIds',
			'ntpUsed': 'ntpIds'}

		tableWriter.writeAttributes(
			aaIds = self.aaIds,
			uncharged_trna_ids = self.uncharged_trna_ids,
			ntpIds = self.ntpIds,
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			gtpPoolSize = self.gtpPoolSize,
			gtpRequestSize = self.gtpRequestSize,
			gtpAllocated = self.gtpAllocated,
			gtpUsed = self.gtpUsed,
			activeRibosomeAllocated = self.activeRibosomeAllocated,
			aaPoolSize = self.aaPoolSize,
			aaRequestSize = self.aaRequestSize,
			aaAllocated = self.aaAllocated,
			aasUsed = self.aasUsed,
			fraction_trna_charged = self.fraction_trna_charged,
			net_charged = self.net_charged,
			ntpPoolSize = self.ntpPoolSize,
			ntpRequestSize = self.ntpRequestSize,
			ntpAllocated = self.ntpAllocated,
			ntpUsed = self.ntpUsed,
			rela_syn = self.rela_syn,
			spot_syn = self.spot_syn,
			spot_deg = self.spot_deg,
			)
