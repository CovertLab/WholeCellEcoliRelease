"""
TranscriptElongationListener
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class TranscriptElongationListener(wholecell.listeners.listener.Listener):
	""" TranscriptElongationListener """

	_name = 'TranscriptElongationListener'

	def __init__(self, *args, **kwargs):
		super(TranscriptElongationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(TranscriptElongationListener, self).initialize(sim, sim_data)

		# Attributes
		self.rnaIds = sim_data.process.transcription.rna_data['id']
		self.attenuated_rnas = self.rnaIds[sim_data.process.transcription.attenuated_rna_indices]
		n_attenuated = len(self.attenuated_rnas)

		# Columns
		self.countRnaSynthesized = np.zeros(sim_data.process.transcription.rna_data.fullArray().size, np.int64)
		self.countNTPsUSed = 0
		self.attenuation_probability = np.zeros(n_attenuated)
		self.counts_attenuated = np.zeros(n_attenuated, np.int64)

	def tableCreate(self, tableWriter):
		subcolumns = {
			'countRnaSynthesized': 'rnaIds',
			'attenuation_probability': 'attenuated_rnas',
			'counts_attenuated': 'attenuated_rnas',
		}

		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaSynthesized = self.countUnits,
			countNTPsUSed = self.countUnits,
			rnaIds = list(self.rnaIds),
			attenuated_rnas = list(self.attenuated_rnas),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countRnaSynthesized = self.countRnaSynthesized,
			countNTPsUSed = self.countNTPsUSed,
			attenuation_probability = self.attenuation_probability,
			counts_attenuated = self.counts_attenuated,
			)
