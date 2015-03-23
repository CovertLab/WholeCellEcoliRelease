#!/usr/bin/env python

"""
ConcentrationChange

@author: John Mason
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 7/18/2014
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class ConcentrationChange(wholecell.listeners.listener.Listener):
	""" ConcentrationChange """

	_name = "ConcentrationChange"

	# Constructor
	def __init__(self, *args, **kwargs):
		super(ConcentrationChange, self).__init__(*args, **kwargs)


	# Construct object graph
	def initialize(self, sim, kb):
		super(ConcentrationChange, self).initialize(sim, kb)

		self.moleculeIDs = kb.process.metabolism.metabolitePoolIDs


	# Allocate memory
	def allocate(self):
		super(ConcentrationChange, self).allocate()

		self.concentrationChange = np.zeros(len(self.moleculeIDs), np.float64)


	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			moleculeIDs = self.moleculeIDs
			)


	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			timeStep = self.timeStep(),
			concentrationChange = self.concentrationChange,
			)
