"""
RnaDegradationListener
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class RnaDegradationListener(wholecell.listeners.listener.Listener):
	""" RnaDegradationListener """

	_name = 'RnaDegradationListener'

	def __init__(self, *args, **kwargs):
		super(RnaDegradationListener, self).__init__(*args, **kwargs)

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaDegradationListener, self).initialize(sim, sim_data)
		self.rnaIds = sim_data.process.transcription.rna_data['id']
		self.cistron_ids = sim_data.process.transcription.cistron_data["id"]
		self.cistron_tu_mapping_matrix = sim_data.process.transcription.cistron_tu_mapping_matrix

	def allocate(self):
		super(RnaDegradationListener, self).allocate()

		self.countRnaDegraded = np.zeros(len(self.rnaIds), np.int64)
		self.count_RNA_degraded_per_cistron = np.zeros(len(self.cistron_ids), np.int64)
		self.nucleotidesFromDegradation = 0
		self.FractionActiveEndoRNases = 0.
		self.DiffRelativeFirstOrderDecay = 0.
		self.FractEndoRRnaCounts = 0.
		self.fragmentBasesDigested = 0

	def update(self):
		self.count_RNA_degraded_per_cistron = self.cistron_tu_mapping_matrix.dot(
			self.countRnaDegraded)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'countRnaDegraded': 'rnaIds',
			'count_RNA_degraded_per_cistron': 'cistron_ids'}

		tableWriter.writeAttributes(
			rnaIds = list(self.rnaIds),
			cistron_ids = list(self.cistron_ids),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countRnaDegraded = self.countRnaDegraded,
			count_RNA_degraded_per_cistron = self.count_RNA_degraded_per_cistron,
			nucleotidesFromDegradation = self.nucleotidesFromDegradation,
			FractionActiveEndoRNases = self.FractionActiveEndoRNases,
			DiffRelativeFirstOrderDecay = self.DiffRelativeFirstOrderDecay,
			FractEndoRRnaCounts = self.FractEndoRRnaCounts,
			fragmentBasesDigested = self.fragmentBasesDigested,
			)
