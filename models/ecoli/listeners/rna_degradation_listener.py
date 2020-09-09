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

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(RnaDegradationListener, self).initialize(sim, sim_data)

		self.countRnaDegraded = np.zeros(sim_data.process.transcription.rna_data.fullArray().size, np.int64)
		self.nucleotidesFromDegradation = 0
		self.FractionActiveEndoRNases = 0.
		self.DiffRelativeFirstOrderDecay = 0.
		self.FractEndoRRnaCounts = 0.
		self.fragmentBasesDigested = 0
		self.rnaIds = sim_data.process.transcription.rna_data['id']

	def tableCreate(self, tableWriter):
		subcolumns = {
			'countRnaDegraded': 'rnaIds'}

		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaDegraded = self.countUnits,
			nucleotidesFromDegradation = self.countUnits,
			FractionActiveEndoRNases = self.countUnits,
			DiffRelativeFirstOrderDecay = self.countUnits,
			FractEndoRRnaCounts = self.countUnits,
			fragmentBasesDigested = self.countUnits,
			rnaIds = list(self.rnaIds),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countRnaDegraded = self.countRnaDegraded,
			nucleotidesFromDegradation = self.nucleotidesFromDegradation,
			FractionActiveEndoRNases = self.FractionActiveEndoRNases,
			DiffRelativeFirstOrderDecay = self.DiffRelativeFirstOrderDecay,
			FractEndoRRnaCounts = self.FractEndoRRnaCounts,
			fragmentBasesDigested = self.fragmentBasesDigested,
			)
