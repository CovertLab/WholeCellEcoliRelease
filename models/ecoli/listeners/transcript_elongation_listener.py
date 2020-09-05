"""
TranscriptElongationListener

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/15
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

		self.countRnaSynthesized = np.zeros(sim_data.process.transcription.rna_data.fullArray().size, np.int64)
		self.countNTPsUSed = 0
		self.rnaIds = sim_data.process.transcription.rna_data['id']


	def tableCreate(self, tableWriter):
		subcolumns = {
			'countRnaSynthesized': 'rnaIds'}

		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaSynthesized = self.countUnits,
			countNTPsUSed = self.countUnits,
			rnaIds = list(self.rnaIds),
			subcolumns = subcolumns)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countRnaSynthesized = self.countRnaSynthesized,
			countNTPsUSed = self.countNTPsUSed,
			)
