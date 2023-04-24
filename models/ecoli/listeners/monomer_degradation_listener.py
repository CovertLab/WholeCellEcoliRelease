"""
MonomerDegradationListener
"""

from __future__ import absolute_import, division, print_function

import numpy as np

import wholecell.listeners.listener

class MonomerDegradationListener(wholecell.listeners.listener.Listener):
	""" MonomerDegradationListener """

	_name = 'MonomerDegradationListener'

	def __init__(self, *args, **kwargs):
		super(MonomerDegradationListener, self).__init__(*args, **kwargs)

		self.countUnits = "counts"

	# Construct object graph
	def initialize(self, sim, sim_data):
		super(MonomerDegradationListener, self).initialize(sim, sim_data)

		self.monomer_ids = sim_data.process.translation.monomer_data['id'].tolist()
		self.monomers_degraded = np.zeros(len(self.monomer_ids), np.int64)

	def tableCreate(self, tableWriter):
		tableWriter.writeAttributes(
			monomer_ids = self.monomer_ids,
			)

	def tableAppend(self, tableWriter):
		tableWriter.append(
			monomers_degraded = self.monomers_degraded,
			)
