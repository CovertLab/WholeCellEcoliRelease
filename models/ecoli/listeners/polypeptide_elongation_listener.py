#!/usr/bin/env python

"""
PolypeptideElongationListener

@author: Derek Macklin
@organization: Covert Lab, Department of Bioengineering, Stanford University
@date: Created 2/1/17
"""

from __future__ import division

import numpy as np

import wholecell.listeners.listener

class PolypeptideElongationListener(wholecell.listeners.listener.Listener):
	""" PolypeptideElongationListener """

	_name = 'PolypeptideElongationListener'

	def __init__(self, *args, **kwargs):
		super(PolypeptideElongationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(PolypeptideElongationListener, self).initialize(sim, sim_data)

		self.countMonomerSynthesized = np.zeros(sim_data.process.translation.monomerData.fullArray().size, np.int64)

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			countMonomerSynthesized = self.countUnits,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			time = self.time(),
			simulationStep = self.simulationStep(),
			countMonomerSynthesized = self.countMonomerSynthesized,
			)
