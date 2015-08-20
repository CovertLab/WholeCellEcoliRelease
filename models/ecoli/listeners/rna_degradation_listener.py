#!/usr/bin/env python

"""
RnaDegradationListener

@author: Javier Carrera
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 1/15/15
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class RnaDegradationListener(wholecell.listeners.listener.Listener):
	""" RnaDegradationListener """

	_name = 'RnaDegradationListener'

	def __init__(self, *args, **kwargs):
		super(RnaDegradationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, kb):
		super(RnaDegradationListener, self).initialize(sim, kb)

		self.countRnaDegraded = np.zeros(kb.process.transcription.rnaData.fullArray().size, np.int64)
		self.nucleotidesFromDegradation = 0
		self.FractionActiveEndoRNases = 0.

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes( # TODO: reconsider attribute names
			countRnaDegraded = self.countUnits,
			nucleotidesFromDegradation = self.countUnits,
			FractionActiveEndoRNases = self.countUnits,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			countRnaDegraded = self.countRnaDegraded,
			nucleotidesFromDegradation = self.nucleotidesFromDegradation,
			FractionActiveEndoRNases = self.FractionActiveEndoRNases,
			)
