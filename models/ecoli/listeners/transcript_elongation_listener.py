#!/usr/bin/env python

"""
TranscriptElongationListener

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/15
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class TranscriptElongationListener(wholecell.listeners.listener.Listener):
	""" TranscriptElongationListener """

	_name = 'TranscriptElongationListener'

	def __init__(self, *args, **kwargs):
		super(TranscriptElongationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(TranscriptElongationListener, self).initialize(sim, kb)

		self.countRnaSynthesized = np.zeros(kb.process.transcription.rnaData.fullArray().size, np.int64)

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaSynthesized = self.countUnits,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			countRnaSynthesized = self.countRnaSynthesized,
			)